"""Fixed rooms used to lock the simulator output."""

from pathlib import Path

import numpy as np

import roomSimSingle
from src.microphone import Microphone
from src.room import Room
from src.room_sim import RoomSim

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
MICRO_CONFIG = ROOT / "micro_config"


def omni_rir():
    room = Room([3.5, 4.0, 2.5], abs_coeff=0.35)
    mic = Microphone(
        [1.2, 1.4, 1.1],
        1,
        orientation=[15.0, -10.0, 5.0],
        direction="omnidirectional",
        micro_config_path=str(MICRO_CONFIG),
    )
    sim = RoomSim(16000, room, [mic], RT60=0.02)
    return sim.create_rir([2.2, 1.7, 1.2], source_off=[30.0, 10.0, 0.0])


def directional_rir():
    frequencies = [125, 250, 500, 1000, 2000, 4000, 8000]
    base = np.array([0.70, 0.75, 0.80, 0.85, 0.90, 0.92, 0.95])
    absorption = [
        base,
        base * 0.98,
        base * 0.95,
        base * 1.02,
        base * 0.90,
        base,
    ]
    mics = [
        Microphone(
            [1.0, 1.2, 1.0],
            1,
            orientation=[20.0, 0.0, 0.0],
            direction="mic_pattern",
            micro_config_path=str(FIXTURES),
        ),
        Microphone(
            [2.5, 3.0, 1.6],
            2,
            orientation=[-30.0, 15.0, 10.0],
            direction="mic_pattern",
            micro_config_path=str(FIXTURES),
        ),
    ]
    room = Room([3.5, 4.0, 2.5], F_abs=frequencies, abs_coeff=absorption)
    sim = RoomSim(16000, room, mics, RT60=None)
    return sim.create_rir(
        [2.0, 2.2, 1.3],
        source_off=[40.0, -20.0, 15.0],
        source_dir=str(FIXTURES / "source_pattern.txt"),
    )


def facade_rir():
    return roomSimSingle.do_everything(
        [4.0, 5.0, 3.0],
        [[1.5, 2.0, 1.2], [2.2, 2.5, 1.4]],
        [3.0, 1.5, 1.6],
        0.03,
    )


def write_directivity_patterns():
    FIXTURES.mkdir(parents=True, exist_ok=True)
    elevation = np.linspace(-90, 90, 181)[:, None] * np.pi / 180
    azimuth = np.linspace(-180, 180, 361)[None, :] * np.pi / 180
    mic = 0.2 + 0.8 * np.cos(elevation) ** 2 * (0.5 + 0.5 * np.cos(azimuth))
    source = 0.3 + 0.7 * np.cos(elevation + 0.4) ** 2 * (0.5 + 0.5 * np.cos(azimuth * 2))
    np.savetxt(FIXTURES / "mic_pattern.txt", mic, fmt="%.18e")
    np.savetxt(FIXTURES / "source_pattern.txt", source, fmt="%.18e")
