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
from matplotlib.patches import FancyBboxPatch
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

INK = "#1c1917"
MUTED = "#78716c"
PAPER = "#f4f1eb"
CARD = "#fffcf8"
TEAL = "#0f766e"
AMBER = "#b45309"
BLUE = "#1d4ed8"
RED = "#9f1239"
LINE = "#e7e5e4"


def _font():
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["axes.unicode_minus"] = False


def _save(fig, path):
    fig.savefig(path, dpi=140, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def _card(fig, x, y, w, h, eyebrow, value, detail, tone):
    ax = fig.add_axes([x, y, w, h])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch(
        (0.0, 0.0), 1.0, 1.0,
        boxstyle="round,pad=0.012,rounding_size=0.04",
        facecolor=CARD,
        edgecolor=LINE,
        linewidth=1.2,
        mutation_aspect=0.55,
        transform=ax.transAxes,
        clip_on=False,
    ))
    ax.add_patch(plt.Rectangle((0.0, 0.0), 0.018, 1.0, color=tone, transform=ax.transAxes, clip_on=False))
    ax.text(0.07, 0.74, eyebrow, color=MUTED, fontsize=10, ha="left", va="center")
    ax.text(0.07, 0.42, value, color=INK, fontsize=18, fontweight="bold", ha="left", va="center")
    ax.text(0.07, 0.16, detail, color=MUTED, fontsize=9, ha="left", va="center")


def draw_scorecard(path, metrics):
    fig = plt.figure(figsize=(13.2, 7.4))
    fig.patch.set_facecolor(PAPER)
    fig.text(0.045, 0.94, "PY_ROOMSIM", color=TEAL, fontsize=12, fontweight="bold", ha="left")
    fig.text(0.045, 0.875, "Оценка качества симулятора", color=INK, fontsize=22, fontweight="bold", ha="left")
    fig.text(
        0.045, 0.825,
        "эталон ≤ 1e-12   ·   прямой звук ±1 отсчёт   ·   RT60 по наклону T20",
        color=MUTED, fontsize=11, ha="left",
    )
    fig.text(
        0.955, 0.94,
        "NumPy %s   ·   SciPy %s   ·   Python %s" % (
            metrics["versions"]["numpy"], metrics["versions"]["scipy"], metrics["versions"]["python"],
        ),
        color=MUTED, fontsize=9, ha="right",
    )

    fixtures = metrics["fixtures"]
    worst = max(row["max_abs"] for row in fixtures)
    omni_path = metrics["direct_path"]
    rt60 = metrics["rt60"]
    repeated = metrics["repeat_identical"]
    rt60_tone = TEAL if rt60["relative_error"] <= 0.35 else AMBER

    cards = [
        ("Эталон", "%.2e" % worst, "макс. |Δ| по трём комнатам", TEAL if all(row["within_atol"] for row in fixtures) else AMBER),
        ("Повторный вызов", "идентичен" if repeated else "различается", "omni, два прогона подряд", TEAL if repeated else AMBER),
        ("Прямой звук", "%+.2f отсч." % omni_path["error_samples"], "пик относительно distance/343", TEAL if omni_path["within_tolerance"] else AMBER),
        ("RT60", "%.3f с" % rt60["estimated_s"], "цель %.2f с, ошибка %.0f%%" % (rt60["target_s"], 100 * rt60["relative_error"]), rt60_tone),
    ]
    gap = 0.018
    width = (0.91 - 3 * gap) / 4
    for index, (eyebrow, value, detail, tone) in enumerate(cards):
        _card(fig, 0.045 + index * (width + gap), 0.52, width, 0.26, eyebrow, value, detail, tone)

    ax = fig.add_axes([0.18, 0.12, 0.76, 0.32])
    ax.set_facecolor(PAPER)
    names = [row["title"] for row in metrics["benchmark"]]
    values = [row["median_ms"] for row in metrics["benchmark"]]
    positions = np.arange(len(names))
    bars = ax.barh(positions, values, color=["#0f766e", "#115e59", "#134e4a"], height=0.62)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=11, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel("медиана create_rir, мс", color=MUTED)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(LINE)
    ax.grid(axis="x", color=LINE)
    ax.set_axisbelow(True)
    for bar, row in zip(bars, metrics["benchmark"]):
        ax.text(
            bar.get_width() + max(values) * 0.02,
            bar.get_y() + bar.get_height() / 2,
            "%.1f мс  ·  %d отсч." % (row["median_ms"], row["samples"]),
            va="center", ha="left", color=INK, fontsize=10,
        )
    ax.set_xlim(0, max(values) * 1.55)
    fig.text(0.18, 0.045, "c = 343 м/с", color=MUTED, fontsize=9)
    _save(fig, path)


def draw_waveforms(path, fixtures):
    fig, axes = plt.subplots(3, 1, figsize=(12.4, 8.2), sharex=False)
    fig.patch.set_facecolor(PAPER)
    fig.suptitle("Импульсные отклики и геометрическая задержка", color=INK, fontsize=16, fontweight="bold")
    colors = (TEAL, AMBER)
    for ax, row in zip(axes, fixtures):
        ax.set_facecolor(CARD)
        actual = row["actual"]
        time_ms = np.arange(actual.shape[0]) / row["sample_rate"] * 1000.0
        for channel in range(actual.shape[1]):
            ax.plot(time_ms, actual[:, channel], color=colors[channel % 2], lw=1.0, label="канал %d" % (channel + 1))
            delay_ms = row["direct_paths"][channel]["expected_samples"] / row["sample_rate"] * 1000.0
            ax.axvline(delay_ms, color=BLUE, ls="--", lw=1.0, alpha=0.9)
        ax.set_ylabel("амплитуда", color=MUTED)
        ax.set_title(row["title"], loc="left", color=INK, fontsize=12)
        ax.tick_params(colors=MUTED)
        for spine in ax.spines.values():
            spine.set_color(LINE)
        ax.grid(color=LINE, lw=0.6)
        ax.legend(frameon=False, fontsize=9)
    axes[-1].set_xlabel("время, мс", color=MUTED)
    fig.text(0.01, 0.01, "Пунктир — задержка distance / 343.", color=MUTED, fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    _save(fig, path)


def draw_fixture_error(path, fixtures):
    fig, axes = plt.subplots(3, 1, figsize=(12.4, 8.0))
    fig.patch.set_facecolor(PAPER)
    fig.suptitle("Модуль разности с эталонными RIR", color=INK, fontsize=16, fontweight="bold")
    for ax, row in zip(axes, fixtures):
        ax.set_facecolor(CARD)
        diff = np.abs(row["actual"] - row["expected"])
        if diff.ndim == 2:
            diff = np.max(diff, axis=1)
        shown = diff.astype(float)
        shown[shown == 0.0] = np.nan
        time_ms = np.arange(diff.size) / row["sample_rate"] * 1000.0
        ax.semilogy(time_ms, shown, color=RED, lw=1.0)
        ax.axhline(FIXTURE_ATOL, color=TEAL, ls="--", lw=1.1, label="порог 1e-12")
        ax.set_ylim(1e-18, 1e-11)
        ax.set_title("%s   max |Δ| = %.2e" % (row["title"], row["max_abs"]), loc="left", color=INK, fontsize=12)
        ax.set_ylabel("|Δ|", color=MUTED)
        ax.tick_params(colors=MUTED)
        for spine in ax.spines.values():
            spine.set_color(LINE)
        ax.grid(color=LINE, which="both", lw=0.5)
        ax.legend(frameon=False, fontsize=9)
    axes[-1].set_xlabel("время, мс", color=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, path)


def draw_energy(path, rt60):
    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    fig.patch.set_facecolor(PAPER)
    ax.set_facecolor(CARD)
    time_s = rt60["time_s"]
    decay = rt60["decay_db"]
    ax.plot(time_s, decay, color=TEAL, lw=1.6, label="кривая Шрёдера")
    mask = rt60["fit_mask"]
    slope = rt60["slope_db_per_s"]
    intercept = rt60["intercept_db"]
    t_line = np.linspace(float(time_s[mask][0]), max(rt60["estimated_s"], float(time_s[mask][-1])), 40)
    ax.plot(t_line, slope * t_line + intercept, color=AMBER, lw=1.4, ls="--", label="экстраполяция T20 к −60 дБ")
    ax.axvline(rt60["target_s"], color=BLUE, ls=":", lw=1.3, label="цель %.2f с" % rt60["target_s"])
    ax.axvline(rt60["estimated_s"], color=AMBER, ls=":", lw=1.2, label="оценка %.3f с" % rt60["estimated_s"])
    ax.set_ylim(-70, 2)
    ax.set_xlim(0, max(float(time_s[-1]), rt60["target_s"], rt60["estimated_s"]) * 1.05)
    ax.set_title("Спад энергии в демо-комнате, RT60 = %.2f с" % rt60["target_s"], loc="left", color=INK, fontsize=15, fontweight="bold")
    ax.set_xlabel("время, с", color=MUTED)
    ax.set_ylabel("уровень, дБ", color=MUTED)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(LINE)
    ax.grid(color=LINE, lw=0.6)
    ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    _save(fig, path)


def _style_axes(ax):
    ax.set_facecolor(CARD)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_color(LINE)


def draw_speech(path, dry, wet, sample_rate):
    fig = plt.figure(figsize=(12.8, 7.8))
    fig.patch.set_facecolor(PAPER)
    grid = fig.add_gridspec(2, 3, width_ratios=(1, 1, 0.045), wspace=0.32, hspace=0.42)
    wave_axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    spec_axes = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    color_ax = fig.add_subplot(grid[1, 2])
    fig.suptitle("Речь до и после свёртки с импульсным откликом", color=INK, fontsize=16, fontweight="bold")
    time_s = np.arange(dry.shape[0]) / sample_rate
    duration = float(time_s[-1]) if time_s.size else 0.0
    peak = max(float(np.max(np.abs(dry))), float(np.max(np.abs(wet))))
    spectrograms = []
    for ax, signal_in, title, color in (
        (wave_axes[0], dry, "Сухой сигнал", TEAL),
        (wave_axes[1], wet, "После свёртки", AMBER),
    ):
        _style_axes(ax)
        ax.plot(time_s, signal_in, color=color, lw=0.6)
        ax.set_xlim(0, duration)
        ax.set_ylim(-peak * 1.05, peak * 1.05)
        ax.set_title(title, loc="left", color=INK)
        ax.set_xlabel("время, с", color=MUTED)
        ax.set_ylabel("амплитуда", color=MUTED)
        ax.grid(color=LINE, lw=0.5)
        frequencies, times, power = signal.spectrogram(signal_in, sample_rate, nperseg=512, noverlap=384)
        spectrograms.append((frequencies, times, 10.0 * np.log10(power + 1e-12)))
    floor = min(float(np.percentile(item[2], 5)) for item in spectrograms)
    ceiling = max(float(np.percentile(item[2], 99.5)) for item in spectrograms)
    image = None
    for ax, (frequencies, times, decibels), title in zip(
        spec_axes, spectrograms, ("Спектр сухого", "Спектр с реверберацией"),
    ):
        _style_axes(ax)
        image = ax.pcolormesh(times, frequencies, decibels, shading="auto", cmap="magma", vmin=floor, vmax=ceiling)
        ax.set_xlim(0, duration)
        ax.set_ylim(0, 8000)
        ax.set_title(title, loc="left", color=INK)
        ax.set_xlabel("время, с", color=MUTED)
        ax.set_ylabel("частота, Гц", color=MUTED)
    fig.colorbar(image, cax=color_ax, label="дБ")
    fig.subplots_adjust(left=0.07, right=0.96, top=0.90, bottom=0.08)
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
