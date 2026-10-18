#!/usr/bin/env python3
"""Build the quality infographic, listening examples, and QUALITY.md."""

import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import soundfile as sf
from scipy import signal

from src.overlap_add_filter import OverlapAddFilter
from tests.quality_metrics import (
    DEMO_MICS,
    DEMO_ROOM_M,
    DEMO_RT60_S,
    DEMO_SOURCE,
    DEMO_SR,
    FIXTURE_ATOL,
    benchmark_scenarios,
    direct_path_report,
    estimate_rt60,
    evaluate_fixtures,
    repeat_is_identical,
    simulate_demo,
)
from utils.audio_loader import AudioLoader

IMG_DIR = ROOT / "img" / "quality"
EXAMPLES_DIR = ROOT / "examples"
LISTENING_PEAK = 0.98

# Okabe–Ito colors, readable in print and for deuteranopia.
BLUE = "#0072B2"
VERMILION = "#D55E00"
BLACK = "#222222"

ROOM_LABELS = {
    "omni": "Omnidirectional",
    "directional": "Directional",
    "facade": "Facade",
}


def _font():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "axes.linewidth": 0.6,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _room_label(row):
    return ROOM_LABELS.get(row["key"], row["title"])


def _panel_title(ax, text):
    ax.set_title(text, loc="left", pad=3, color="black")


def draw_scorecard(path, metrics):
    fixtures = metrics["fixtures"]
    benchmark = metrics["benchmark"]
    direct = metrics["direct_path"]
    rt60 = metrics["rt60"]
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 3.05))

    error_ax, time_ax = axes
    labels = [_room_label(row) for row in fixtures]
    positions = np.arange(len(labels))
    errors = [row["max_abs"] for row in fixtures]
    error_ax.barh(positions, errors, color=BLUE, height=0.55, zorder=2)
    error_ax.set_xscale("log")
    error_ax.set_xlim(1e-16, 3e-11)
    error_ax.axvline(FIXTURE_ATOL, color=BLACK, ls="--", lw=0.7, zorder=1)
    error_ax.text(
        FIXTURE_ATOL, 1.03,
        r"$10^{-12}$",
        transform=error_ax.get_xaxis_transform(),
        ha="center", va="bottom", fontsize=7, clip_on=False,
    )
    error_ax.set_yticks(positions)
    error_ax.set_yticklabels(labels)
    error_ax.invert_yaxis()
    error_ax.set_xlabel(r"Maximum $|\Delta|$")
    _panel_title(error_ax, "(a)  Deviation from the reference")
    for y_pos, value in zip(positions, errors):
        error_ax.text(value * 1.8, y_pos, "%.1e" % value, va="center", ha="left", fontsize=7, color=BLACK)

    times = [row["median_ms"] for row in benchmark]
    time_ax.barh(positions, times, color=BLUE, height=0.55)
    time_ax.set_yticks(positions)
    time_ax.set_yticklabels(labels)
    time_ax.invert_yaxis()
    time_ax.set_xlabel(r"Median runtime (ms)")
    time_ax.set_xlim(0, max(times) * 1.38)
    _panel_title(time_ax, "(b)  Runtime of create_rir")
    for y_pos, row in zip(positions, benchmark):
        time_ax.text(
            row["median_ms"] + max(times) * 0.03,
            y_pos,
            "%.1f ms" % row["median_ms"],
            va="center", ha="left", fontsize=7, color=BLACK,
        )

    repeat = "identical" if metrics["repeat_identical"] else "not identical"
    fig.text(
        0.0, 0.02,
        "Omnidirectional direct-path error: %+.2f samples.\n"
        r"Schroeder $T_{20}$ estimate: %.3f s (target %.2f s, relative error %.0f%%). Repeated call: %s."
        % (direct["error_samples"], rt60["estimated_s"], rt60["target_s"], 100 * rt60["relative_error"], repeat),
        fontsize=7, ha="left", va="bottom", color=BLACK, linespacing=1.35,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    _save(fig, path)


def draw_waveforms(path, fixtures):
    fig, axes = plt.subplots(3, 1, figsize=(7.16, 6.15), sharex=False)
    letters = ("(a)", "(b)", "(c)")
    for ax, row, letter in zip(axes, fixtures, letters):
        actual = row["actual"]
        time_ms = np.arange(actual.shape[0]) / row["sample_rate"] * 1000.0
        delay_labeled = False
        for channel in range(actual.shape[1]):
            label = "Impulse response" if actual.shape[1] == 1 else "Channel %d" % (channel + 1)
            ax.plot(time_ms, actual[:, channel], color=(BLUE, VERMILION)[channel % 2], lw=0.7, label=label)
            delay_ms = row["direct_paths"][channel]["expected_samples"] / row["sample_rate"] * 1000.0
            ax.axvline(
                delay_ms, color="0.25", ls="--", lw=0.6,
                label="Geometric delay" if not delay_labeled else None,
            )
            delay_labeled = True
        ax.set_ylabel("Amplitude")
        _panel_title(ax, "%s  %s" % (letter, _room_label(row)))
        ax.legend(loc="upper right", handlelength=1.6)
    axes[-1].set_xlabel("Time (ms)")
    fig.tight_layout()
    _save(fig, path)


def draw_fixture_error(path, fixtures):
    fig, axes = plt.subplots(3, 1, figsize=(7.16, 6.15))
    letters = ("(a)", "(b)", "(c)")
    for index, (ax, row, letter) in enumerate(zip(axes, fixtures, letters)):
        diff = np.abs(row["actual"] - row["expected"])
        if diff.ndim == 2:
            diff = np.max(diff, axis=1)
        shown = diff.astype(float)
        shown[shown == 0.0] = np.nan
        time_ms = np.arange(diff.size) / row["sample_rate"] * 1000.0
        ax.semilogy(time_ms, shown, color=BLACK, lw=0.7)
        ax.axhline(FIXTURE_ATOL, color=VERMILION, ls="--", lw=0.7, label=r"Tolerance $10^{-12}$")
        ax.set_ylim(1e-18, 1e-11)
        ax.set_ylabel(r"$|\Delta|$")
        _panel_title(ax, r"%s  %s,  $\max|\Delta| = %.2e$" % (letter, _room_label(row), row["max_abs"]))
        if index == 0:
            ax.legend(loc="upper right")
    axes[-1].set_xlabel("Time (ms)")
    fig.tight_layout()
    _save(fig, path)


def draw_energy(path, rt60):
    fig, ax = plt.subplots(figsize=(7.16, 3.35))
    time_s = rt60["time_s"]
    decay = rt60["decay_db"]
    ax.plot(time_s, decay, color=BLACK, lw=1.0, label="Schroeder curve")
    mask = rt60["fit_mask"]
    slope = rt60["slope_db_per_s"]
    intercept = rt60["intercept_db"]
    t_line = np.linspace(float(time_s[mask][0]), max(rt60["estimated_s"], float(time_s[mask][-1])), 50)
    ax.plot(t_line, slope * t_line + intercept, color=VERMILION, lw=0.9, ls="--", label=r"$T_{20}$ extrapolated to $-60$ dB")
    ax.axvline(rt60["target_s"], color=BLUE, ls=":", lw=0.8, label="Target %.2f s" % rt60["target_s"])
    ax.axvline(rt60["estimated_s"], color=VERMILION, ls=":", lw=0.8, label="Estimate %.3f s" % rt60["estimated_s"])
    ax.set_ylim(-70, 2)
    ax.set_xlim(0, max(float(time_s[-1]), rt60["target_s"], rt60["estimated_s"]) * 1.02)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Level (dB)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.38, 0.98), frameon=False, handlelength=1.8, borderaxespad=0.0)
    fig.tight_layout()
    _save(fig, path)


def draw_speech(path, dry, wet, sample_rate):
    fig = plt.figure(figsize=(7.16, 5.35))
    grid = fig.add_gridspec(2, 3, width_ratios=(1, 1, 0.045), wspace=0.38, hspace=0.48)
    wave_axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    spec_axes = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    color_ax = fig.add_subplot(grid[1, 2])
    time_s = np.arange(dry.shape[0]) / sample_rate
    duration = float(time_s[-1]) if time_s.size else 0.0
    peak = max(float(np.max(np.abs(dry))), float(np.max(np.abs(wet))))
    spectrograms = []
    for ax, signal_in, title, color in (
        (wave_axes[0], dry, "(a)  Dry signal", BLACK),
        (wave_axes[1], wet, "(b)  Reverberant signal", VERMILION),
    ):
        ax.plot(time_s, signal_in, color=color, lw=0.35)
        ax.set_xlim(0, duration)
        ax.set_ylim(-peak * 1.05, peak * 1.05)
        _panel_title(ax, title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        frequencies, times, power = signal.spectrogram(signal_in, sample_rate, nperseg=512, noverlap=384)
        spectrograms.append((frequencies, times, 10.0 * np.log10(power + 1e-12)))
    floor = min(float(np.percentile(item[2], 5)) for item in spectrograms)
    ceiling = max(float(np.percentile(item[2], 99.5)) for item in spectrograms)
    image = None
    for ax, (frequencies, times, decibels), title in zip(
        spec_axes, spectrograms, ("(c)  Dry spectrogram", "(d)  Reverberant spectrogram"),
    ):
        image = ax.pcolormesh(
            times, frequencies, decibels, shading="auto", cmap="cividis", vmin=floor, vmax=ceiling, rasterized=True,
        )
        ax.set_xlim(0, duration)
        ax.set_ylim(0, sample_rate / 2.0)
        _panel_title(ax, title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Frequency (Hz)")
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label("Level (dB)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7, width=0.6, length=3)
    fig.subplots_adjust(left=0.09, right=0.94, top=0.96, bottom=0.08)
    _save(fig, path)


def write_listening_wav(path, data, sample_rate):
    peak = float(np.max(np.abs(data)))
    scaled = np.zeros_like(data) if peak == 0.0 else np.asarray(data, dtype=float) / peak * LISTENING_PEAK
    sf.write(path, scaled, sample_rate, subtype="PCM_16")
    return peak


def render_markdown(metrics):
    fixtures = "\n".join(
        "| %s | %.3e | %s |" % (row["title"], row["max_abs"], "да" if row["within_atol"] else "нет")
        for row in metrics["fixtures"]
    )
    timings = "\n".join(
        "| %s | %.1f | %d |" % (row["title"], row["median_ms"], row["samples"])
        for row in metrics["benchmark"]
    )
    path = metrics["direct_path"]
    rt60 = metrics["rt60"]
    versions = metrics["versions"]
    return """# Оценка качества PY_ROOMSIM

Отчёт собран скриптом `scripts/make_quality_report.py` {generated}. Сравниваются три эталонные комнаты из `tests/fixtures/`, геометрическая задержка прямого звука и оценка RT60 по кривой Шрёдера. Побитовое равенство RIR на этой платформе не выполняется: расхождение сидит в последних битах FFT и `lfilter`, поэтому порог качества — `1e-12`.

Версии прогона: Python {python}, NumPy {numpy}, SciPy {scipy}, soundfile {soundfile}.

## Сводка

![Сводка качества](img/quality/scorecard.png)

| Проверка | Результат |
| --- | --- |
| Максимум \\|Δ\\| по эталонам | {worst:.3e} |
| Повторный вызов omni | {repeat} |
| Ошибка задержки прямого звука | {error:+.2f} отсчёта (пик {peak}, геометрия {expected:.2f}) |
| RT60, цель {target:.2f} с | оценка {estimated:.3f} с, относительная ошибка {rel:.0f}% |
| Поглощение демо-комнаты | {absorption:.4f} |

## Импульсные отклики

Пунктир на графике — задержка `distance / 343`. Для всенаправленного сценария допуск проверки — один отсчёт.

![Импульсные отклики](img/quality/rir_waveforms.png)

## Расхождение с эталоном

| Комната | max \\|Δ\\| | в пределах 1e-12 |
| --- | --- | --- |
{fixtures}

![Ошибка относительно эталона](img/quality/fixture_error.png)

## Спад энергии

Демо-комната {room} м, RT60 {target:.2f} с, микрофоны всенаправленный и кардиоида. Наклон снят на участке Шрёдера от −5 до −25 дБ и экстраполирован к −60 дБ (T20). Это проверка порядка величины: хвост импульса обрезан длиной, которую задаёт RT60, а воздух в симуляторе добавляет затухание сверх формулы Норриса–Эйринга.

![Спад энергии](img/quality/energy_decay.png)

## Пример прослушивания

`wav/4.wav` приведён к 16 кГц через `AudioLoader` и свёрнут с первым каналом демо-импульса функцией `OverlapAddFilter`. Свёртка возвращает сигнал той же длины, что и сухой фрагмент. Файлы ниже нормированы к пику {listening:.2f} для прослушивания; график показывает амплитуду до нормировки.

- [examples/rir_stereo.wav](examples/rir_stereo.wav) — стереоимпульс, сырой пик {rir_peak:.3f}
- [examples/speech_reverb.wav](examples/speech_reverb.wav) — речь после свёртки, сырой пик {wet_peak:.3f}, сухой пик {dry_peak:.3f}

![Речь до и после](img/quality/speech_before_after.png)

## Время расчёта

Медиана трёх вызовов `create_rir` на эталонных сценариях.

| Комната | медиана, мс | отсчёты |
| --- | --- | --- |
{timings}

## Как повторить

```bash
python3 -m pip install --user -r requirements-report.txt
python3 -m unittest tests.test_quality_metrics -v
python3 scripts/make_quality_report.py
```

Числа этого файла и `img/quality/metrics.json` перезаписываются скриптом.
""".format(
        generated=metrics["generated_at"],
        python=versions["python"],
        numpy=versions["numpy"],
        scipy=versions["scipy"],
        soundfile=versions["soundfile"],
        worst=max(row["max_abs"] for row in metrics["fixtures"]),
        repeat="идентичен" if metrics["repeat_identical"] else "различается",
        error=path["error_samples"],
        peak=path["peak_index"],
        expected=path["expected_samples"],
        target=rt60["target_s"],
        estimated=rt60["estimated_s"],
        rel=100 * rt60["relative_error"],
        absorption=rt60["absorption"],
        room="×".join("%.1f" % value for value in DEMO_ROOM_M),
        fixtures=fixtures,
        listening=LISTENING_PEAK,
        rir_peak=metrics["examples"]["rir_stereo"]["raw_peak"],
        wet_peak=metrics["examples"]["speech_reverb"]["raw_peak"],
        dry_peak=metrics["examples"]["speech_reverb"]["dry_peak"],
        timings=timings,
    )


def public_metrics(metrics):
    exported = json.loads(json.dumps(metrics, default=_json_default))
    for row in exported["fixtures"]:
        row.pop("actual", None)
        row.pop("expected", None)
    rt60 = exported["rt60"]
    for key in ("decay_db", "time_s", "fit_mask"):
        rt60.pop(key, None)
    return exported


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value))


def main():
    _font()
    print("fixtures", flush=True)
    fixtures = evaluate_fixtures()
    print("repeat", flush=True)
    repeated = repeat_is_identical()
    print("benchmark", flush=True)
    benchmark = benchmark_scenarios()
    print("demo room", flush=True)
    demo, absorption = simulate_demo()
    rt60_fit = estimate_rt60(demo, DEMO_SR)
    relative_error = abs(rt60_fit["rt60_s"] - DEMO_RT60_S) / DEMO_RT60_S
    print("speech", flush=True)
    sample_rate, dry = AudioLoader().load_audio(str(ROOT / "wav" / "4.wav"))
    wet = OverlapAddFilter(demo[:, 0], dry)

    IMG_DIR.mkdir(parents=True, exist_ok=True)
    EXAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    rir_peak = write_listening_wav(EXAMPLES_DIR / "rir_stereo.wav", demo, DEMO_SR)
    wet_peak = write_listening_wav(EXAMPLES_DIR / "speech_reverb.wav", wet, sample_rate)

    omni = next(row for row in fixtures if row["key"] == "omni")
    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "soundfile": sf.__version__,
    }
    metrics = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "versions": versions,
        "fixture_atol": FIXTURE_ATOL,
        "fixtures": fixtures,
        "repeat_identical": repeated,
        "direct_path": omni["direct_paths"][0],
        "rt60": {
            "target_s": DEMO_RT60_S,
            "estimated_s": rt60_fit["rt60_s"],
            "relative_error": float(relative_error),
            "absorption": absorption,
            "room_m": list(DEMO_ROOM_M),
            "source_m": list(DEMO_SOURCE),
            "mics": [{"pos": list(spec["pos"]), "direction": spec["direction"]} for spec in DEMO_MICS],
            "sample_rate": DEMO_SR,
            "samples": int(demo.shape[0]),
            "slope_db_per_s": rt60_fit["slope_db_per_s"],
            "fit_points": rt60_fit["fit_points"],
            "time_s": rt60_fit["time_s"],
            "decay_db": rt60_fit["decay_db"],
            "fit_mask": rt60_fit["fit_mask"],
            "intercept_db": rt60_fit["intercept_db"],
        },
        "benchmark": benchmark,
        "examples": {
            "rir_stereo": {
                "path": "examples/rir_stereo.wav",
                "sample_rate": DEMO_SR,
                "samples": int(demo.shape[0]),
                "channels": int(demo.shape[1]),
                "raw_peak": rir_peak,
            },
            "speech_reverb": {
                "path": "examples/speech_reverb.wav",
                "sample_rate": int(sample_rate),
                "samples": int(wet.shape[0]),
                "raw_peak": wet_peak,
                "dry_peak": float(np.max(np.abs(dry))),
            },
        },
    }
    # Keep the geometric check for the demo impulse in the json as well.
    metrics["rt60"]["direct_path"] = direct_path_report(demo, DEMO_MICS[0]["pos"], DEMO_SOURCE, DEMO_SR)

    print("plots", flush=True)
    draw_scorecard(IMG_DIR / "scorecard.png", metrics)
    draw_waveforms(IMG_DIR / "rir_waveforms.png", fixtures)
    draw_fixture_error(IMG_DIR / "fixture_error.png", fixtures)
    draw_energy(IMG_DIR / "energy_decay.png", metrics["rt60"])
    draw_speech(IMG_DIR / "speech_before_after.png", dry, wet, sample_rate)

    (IMG_DIR / "metrics.json").write_text(
        json.dumps(public_metrics(metrics), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (ROOT / "QUALITY.md").write_text(render_markdown(public_metrics(metrics)), encoding="utf-8")
    print("wrote", IMG_DIR, "and", EXAMPLES_DIR, flush=True)
    print("rt60", metrics["rt60"]["estimated_s"], "rel", metrics["rt60"]["relative_error"], flush=True)


if __name__ == "__main__":
    main()
