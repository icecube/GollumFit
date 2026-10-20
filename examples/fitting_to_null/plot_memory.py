#!/usr/bin/env python3
"""Summarise the memory study (bench_memory.slurm -> memory_<job>.out, MEMJSON lines).

Prints a table of resident memory per stage for every run, fits memory = a + b * (simulated events)
per build and thread count, and writes plots/memory_scaling.{pdf,png}.
Usage: python plot_memory.py memory_<job>.out
"""
import json, os, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
plt.style.use(os.path.join(here, "figures.mplstyle"))
if subprocess.call(["which", "latex"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
    plt.rcParams.update({"text.usetex": False, "font.family": "serif", "mathtext.fontset": "cm"})
N1 = 334683
runs = []
for line in open(sys.argv[1]):
    if line.startswith("MEMJSON "):
        r = json.loads(line[8:])
        build, kk, tt = r["label"].split()
        r["build"], r["k"], r["t"] = build, int(kk[1:]), int(tt[1:])
        r["n"] = N1 * r["k"]
        r["s"] = {s["stage"]: s for s in r["stages"]}
        runs.append(r)
if not runs:
    sys.exit("no MEMJSON lines")

KEYS = [("constructor (compact file, splines, histogram)", "after load", "rss_mb"),
        ("data + likelihood problem", "+ data, problem", "rss_mb"),
        ("first value evaluation", "after 1st eval", "rss_mb"),
        ("short fit (value + gradient)", "peak (whole run)", "peak_mb")]
print("%-16s %10s" % ("run", "events") + "".join("%18s" % k[1] for k in KEYS))
for r in sorted(runs, key=lambda r: (r["build"], r["t"], r["k"])):
    print("%-16s %10d" % (r["label"], r["n"]) + "".join("%15.0f MB" % r["s"][k][f] for k, _, f in KEYS))

print("\nlinear fits  memory = intercept + slope * events")
fits = {}
for build in ("original", "fast"):
    for t in (1, 8):
        sel = sorted([r for r in runs if r["build"] == build and r["t"] == t], key=lambda r: r["n"])
        if len(sel) < 2:
            continue
        n = np.array([r["n"] for r in sel], float)
        for k, short, f in KEYS:
            y = np.array([r["s"][k][f] for r in sel])
            b, a = np.polyfit(n, y, 1)
            fits[(build, t, short)] = (a, b * 1024 ** 2)
            print("  %-8s %d thread(s)  %-17s intercept %6.0f MB   slope %5.0f bytes/event" % (build, t, short, a, b * 1024 ** 2))

print("\nfast - original, per event (bytes):")
for t in (1, 8):
    for _, short, _ in KEYS:
        if (("fast", t, short) in fits) and (("original", t, short) in fits):
            d = fits[("fast", t, short)][1] - fits[("original", t, short)][1]
            da = fits[("fast", t, short)][0] - fits[("original", t, short)][0]
            print("  %d thread(s) %-17s %+6.0f bytes/event, %+5.0f MB fixed" % (t, short, d, da))

fig, axes = plt.subplots(1, 2, figsize=(8, 3.3), sharey=True)
for ax, (key, short, f) in zip(axes, (KEYS[2], KEYS[3])):
    for build, color in (("original", "#4C72B0"), ("fast", "#DD8452")):
        for t, ls, mk in ((1, "-", "o"), (8, "--", "s")):
            sel = sorted([r for r in runs if r["build"] == build and r["t"] == t], key=lambda r: r["n"])
            if not sel:
                continue
            n = np.array([r["n"] for r in sel]) / 1e6
            y = np.array([r["s"][key][f] for r in sel]) / 1024
            ax.plot(n, y, ls, marker=mk, color=color, ms=4, lw=1.5,
                    label="%s, %d thread%s" % (build, t, "s" if t > 1 else ""), mfc=color if t == 1 else "none")
    ax.set_xlabel("Simulated events [millions]", fontsize=11)
    ax.text(0.97, 0.05, {"after 1st eval": "resident after first evaluation", "peak (whole run)": "peak over the whole run"}[short],
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
axes[0].set_ylabel("Memory [GB]", fontsize=11)
axes[0].legend(fontsize=8, loc="upper left")
for ax in axes:
    ax.tick_params(which="minor", length=3)
fig.tight_layout()
os.makedirs(os.path.join(here, "plots"), exist_ok=True)
fig.savefig(os.path.join(here, "plots", "memory_scaling.pdf")); fig.savefig(os.path.join(here, "plots", "memory_scaling.png"))
print("\nwrote plots/memory_scaling.pdf")
