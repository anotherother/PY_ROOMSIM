"""Quality measures shared by the unit tests and the infographic report.

The omni, directional, and facade geometries match ``tests.scenarios``.
The speed of sound matches ``src.room_sim`` (343 m/s).
"""

import os
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from src.acoustics import rt60_to_absorption
from src.microphone import Microphone
from src.room import Room
from src.room_sim import RoomSim
from tests.scenarios import FIXTURES, directional_rir, facade_rir, omni_rir

ROOT = Path(__file__).resolve().parents[1]
SPEED_OF_SOUND = 343.0
FIXTURE_ATOL = 1e-12
DIRECT_PATH_TOLERANCE_SAMPLES = 1.0

DEMO_ROOM_M = (3.5, 4.0, 2.5)
DEMO_RT60_S = 0.3
DEMO_SR = 16000
DEMO_SOURCE = (2.0, 2.0, 1.3)
DEMO_MICS = (
    {"pos": (1.2, 1.5, 1.2), "direction": "omnidirectional"},
    {"pos": (2.4, 2.6, 1.4), "direction": "cardioid"},
)

SCENARIOS = (
    {
        "key": "omni",
        "title": "Всенаправленный",
        "fixture": "rir_omni.npy",
        "simulate": omni_rir,
        "sample_rate": 16000,
        "channels": (
            {"mic": (1.2, 1.4, 1.1), "source": (2.2, 1.7, 1.2)},
        ),
    },
    {
        "key": "directional",
        "title": "Направленный",
        "fixture": "rir_directional.npy",
        "simulate": directional_rir,
        "sample_rate": 16000,
        "channels": (
            {"mic": (1.0, 1.2, 1.0), "source": (2.0, 2.2, 1.3)},
            {"mic": (2.5, 3.0, 1.6), "source": (2.0, 2.2, 1.3)},
        ),
    },
    {
        "key": "facade",
        "title": "Фасад",
        "fixture": "rir_facade.npy",
        "simulate": facade_rir,
        "sample_rate": 16000,
        "channels": (
            {"mic": (1.5, 2.0, 1.2), "source": (3.0, 1.5, 1.6)},
            {"mic": (2.2, 2.5, 1.4), "source": (3.0, 1.5, 1.6)},
        ),
    },
)


@contextmanager
def at_repo_root():
    previous = os.getcwd()
    os.chdir(ROOT)
    try:
        yield
    finally:
        os.chdir(previous)


def direct_path_report(rir, mic, source, sample_rate, channel=0):
    """Compare the loudest sample with the geometric delay ``distance / 343``."""
    distance = float(np.linalg.norm(np.asarray(source, dtype=float) - np.asarray(mic, dtype=float)))
    delay = distance / SPEED_OF_SOUND * sample_rate
    column = np.asarray(rir)
    if column.ndim == 2:
        column = column[:, channel]
    peak = int(np.argmax(np.abs(column)))
    error = float(peak - delay)
    return {
        "distance_m": distance,
        "expected_samples": float(delay),
        "peak_index": peak,
        "error_samples": error,
        "within_tolerance": abs(error) <= DIRECT_PATH_TOLERANCE_SAMPLES,
        "peak_value": float(column[peak]),
    }


def evaluate_fixtures():
    """Simulate the three reference rooms and compare them with the ``.npy`` fixtures."""
    rows = []
    with at_repo_root():
        for spec in SCENARIOS:
            actual = spec["simulate"]()
            expected = np.load(FIXTURES / spec["fixture"])
            max_abs = float(np.max(np.abs(actual - expected)))
            rows.append({
                "key": spec["key"],
                "title": spec["title"],
                "max_abs": max_abs,
                "within_atol": max_abs <= FIXTURE_ATOL,
                "samples": int(actual.shape[0]),
                "channels": int(actual.shape[1] if actual.ndim == 2 else 1),
                "direct_paths": [
                    direct_path_report(actual, item["mic"], item["source"], spec["sample_rate"], channel)
                    for channel, item in enumerate(spec["channels"])
                ],
                "actual": actual,
                "expected": expected,
                "sample_rate": spec["sample_rate"],
            })
    return rows


def repeat_is_identical():
    with at_repo_root():
        first = omni_rir()
        second = omni_rir()
    return bool(np.array_equal(first, second))


def benchmark_scenarios(repeats=3):
    rows = []
    with at_repo_root():
        for spec in SCENARIOS:
            elapsed_ms = []
            last = None
            for _ in range(repeats):
                started = time.perf_counter()
                last = spec["simulate"]()
                elapsed_ms.append((time.perf_counter() - started) * 1000.0)
            rows.append({
                "key": spec["key"],
                "title": spec["title"],
                "median_ms": float(np.median(elapsed_ms)),
                "min_ms": float(np.min(elapsed_ms)),
                "samples": int(last.shape[0]),
                "channels": int(last.shape[1] if last.ndim == 2 else 1),
            })
    return rows


def schroeder_decay_db(channel):
    energy = np.square(np.asarray(channel, dtype=float))
    decay = np.cumsum(energy[::-1])[::-1]
    reference = float(decay[0])
    if reference <= 0.0:
        raise ValueError("impulse energy is zero")
    return 10.0 * np.log10(np.maximum(decay, np.finfo(float).tiny) / reference)


def estimate_rt60(rir, sample_rate, channel=0, fit_from_db=-5.0, fit_to_db=-25.0):
    """T20 extrapolated to -60 dB on the Schroeder decay of one channel."""
    column = np.asarray(rir)
    if column.ndim == 2:
        column = column[:, channel]
    decay_db = schroeder_decay_db(column)
    time_s = np.arange(decay_db.size, dtype=float) / sample_rate
    mask = (decay_db <= fit_from_db) & (decay_db >= fit_to_db)
    if int(np.count_nonzero(mask)) < 8:
        raise RuntimeError("Schroeder curve does not cover the T20 fit range")
    slope, intercept = np.polyfit(time_s[mask], decay_db[mask], 1)
    if slope >= 0.0:
        raise RuntimeError("Schroeder slope is not a decay")
    return {
        "rt60_s": float(-60.0 / slope),
        "slope_db_per_s": float(slope),
        "intercept_db": float(intercept),
        "fit_from_db": float(fit_from_db),
        "fit_to_db": float(fit_to_db),
        "fit_points": int(np.count_nonzero(mask)),
        "decay_db": decay_db,
        "time_s": time_s,
        "fit_mask": mask,
    }


def simulate_demo():
    """Longer shoebox used for the reverberation-time curve and the listening example."""
    absorption = float(np.asarray(rt60_to_absorption(list(DEMO_ROOM_M), DEMO_RT60_S)))
    room = Room(list(DEMO_ROOM_M), abs_coeff=absorption)
    mics = [
        Microphone(
            list(spec["pos"]),
            index,
            orientation=[0.0, 0.0, 0.0],
            direction=spec["direction"],
            micro_config_path=str(ROOT / "micro_config"),
        )
        for index, spec in enumerate(DEMO_MICS, start=1)
    ]
    simulator = RoomSim(DEMO_SR, room, mics, RT60=DEMO_RT60_S)
    rir = simulator.create_rir(list(DEMO_SOURCE))
    return rir, absorption
