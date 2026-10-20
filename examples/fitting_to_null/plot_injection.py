#!/usr/bin/env python3
"""Summarise the injection-recovery fits (fit_injection.py via bench_injection.slurm).

Reads every inject_*.out in this directory (one or more fit blocks each, starting with the
'build ... priors ... data ...' line written by bench_injection.slurm), the injected values from the
data files, and writes to plots/:
  injection_recovery_<data>.pdf   (fit - truth) / prior width for both builds, one panel per prior setting
  injection_summary.pdf           table: likelihood at truth and at the best fit, evaluations, time, largest deviation
and prints the same numbers. Run with PYTHONNOUSERSITE=1 on fasrc.
"""
import glob, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import bench_common as bc

plt.style.use(os.path.join(here, "figures.mplstyle"))
if subprocess.call(["which", "latex"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
    plt.rcParams.update({"text.usetex": False, "font.family": "serif", "mathtext.fontset": "cm"})
os.makedirs(os.path.join(here, "plots"), exist_ok=True)
C_OLD, C_NEW = "#4C72B0", "#DD8452"
BUILD_LABEL = {"head": "original", "new": "fast"}
DATA_LABEL = {"injection_corr_data": "correlated draw", "injection_data": "independent draw"}


def parse_blocks(path):
    blocks, cur = [], None
    for line in open(path):
        line = line.strip()
        m = re.match(r"^build (\S+) priors (\S+)(?: data (\S+))? host (\S+)", line)
        if m:
            cur = {"build": m.group(1), "priors": m.group(2), "data": m.group(3) or "injection_data.npz",
                   "host": m.group(4), "params": {}, "file": os.path.basename(path)}
            blocks.append(cur)
            continue
        if cur is None:
            continue
        if line.startswith("LLH at truth:"):
            cur["llh_truth"] = float(line.split()[3])
        elif line.startswith("MinLLH wall time:"):
            cur["wall"] = float(line.split()[3])
        elif line.startswith("llh:"):
            cur["llh"] = float(line.split()[1])
        elif line.startswith("nEval:"):
            cur["neval"] = int(line.split()[1])
        elif line.startswith("succeeded:"):
            cur["succeeded"] = line.split()[1] == "True"
        else:
            mm = re.match(r"^([A-Za-z0-9]+) (-?[0-9.eE+-]+)$", line)
            if mm and mm.group(1) in bc.SYST:
                cur["params"][mm.group(1)] = float(mm.group(2))
    return [b for b in blocks if "llh" in b and len(b["params"]) == len(bc.NAMES)]


fits = {}
for f in sorted(glob.glob(os.path.join(here, "inject_*.out")), key=os.path.getmtime):
    for b in parse_blocks(f):
        fits[(b["data"], b["priors"], b["build"])] = b  # latest wins
if not fits:
    sys.exit("no finished injection fits found")

names = bc.NAMES
center = np.array([bc.SYST[n][1] for n in names])
width = np.array([bc.SYST[n][2] for n in names])


def fmt_duration(t):
    if t >= 3600:
        return "%d h %d min" % (t // 3600, round((t % 3600) / 60))
    if t >= 60:
        return "%d min %d s" % (t // 60, round(t % 60))
    return "%.0f s" % t


rows, cmp_rows = [], []
for data in sorted({k[0] for k in fits}):
    d = np.load(os.path.join(here, data))
    truth = np.array([dict(zip([str(n) for n in d["names"]], d["truth"]))[n] for n in names])
    tag = os.path.splitext(data)[0]
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.6), sharex=True, gridspec_kw={"hspace": 0.12})
    x = np.arange(len(names))
    for ax, priors in zip(axes, ("truth", "nominal")):
        if priors == "nominal":
            ax.plot(x, (center - truth) / width, "x", color="0.55", ms=4, label="prior centre (nominal)")
        for build, color, marker, off in (("head", C_OLD, "o", -0.15), ("new", C_NEW, "s", 0.15)):
            b = fits.get((data, priors, build))
            if b is None:
                continue
            p = np.array([b["params"][n] for n in names])
            dev = (p - truth) / width
            ax.plot(x + off, dev, marker, color=color, ms=3.5, ls="none", mfc=color if build == "head" else "none", mew=1.0,
                    label=r"%s: $-\ln L$ %.4f" % (BUILD_LABEL[build], b["llh"]))
            k = int(np.argmax(np.abs(dev)))
            rows.append((tag, priors, BUILD_LABEL[build], b["llh_truth"], b["llh"], b["llh"] - b["llh_truth"],
                         b["neval"], fmt_duration(b["wall"]), "%.3f (%s)" % (abs(dev[k]), names[k]), b.get("succeeded")))
        bo, bn = fits.get((data, priors, "head")), fits.get((data, priors, "new"))
        if bo and bn:
            dp = (np.array([bn["params"][n] for n in names]) - np.array([bo["params"][n] for n in names])) / width
            k = int(np.argmax(np.abs(dp)))
            cmp_rows.append((tag, priors, bn["llh"] - bo["llh"], "%.3f (%s)" % (abs(dp[k]), names[k])))
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_ylabel("(fit $-$ truth) / prior width\npriors centred at %s" % priors, fontsize=10)
        ax.legend(fontsize=7.5, loc="upper left", ncol=3, frameon=False)
        ax.tick_params(axis="x", which="both", length=0)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(names, rotation=90, fontsize=7)
    fig.savefig(os.path.join(here, "plots", "injection_recovery_%s.pdf" % tag))
    fig.savefig(os.path.join(here, "plots", "injection_recovery_%s.png" % tag))
    plt.close(fig)

# ---------------- summary table ----------------
header = ("injection", "priors at", "build", r"$-\ln L$ at truth", r"$-\ln L$ best fit", "best $-$ truth", "evals", "fit time", "max $|$fit$-$truth$|$/width")
cols = [(0.0, "left"), (0.15, "left"), (0.24, "left"), (0.42, "right"), (0.56, "right"), (0.655, "right"), (0.71, "right"), (0.80, "right"), (1.0, "right")]
n = len(rows) + 1.5 + (len(cmp_rows) + 1.5 if cmp_rows else 0)
fig, ax = plt.subplots(figsize=(10.5, 0.27 * n + 0.3))
ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, n)
y = n - 0.75
for (xx, ha), t in zip(cols, header):
    ax.text(xx, y, t, ha=ha, va="center", fontsize=8.5)
ax.plot([0, 1], [n, n], color="black", lw=1.2); ax.plot([0, 1], [n - 1.5, n - 1.5], color="black", lw=0.6)
y = n - 2.0
for r in rows:
    vals = (DATA_LABEL.get(r[0], r[0].replace("_", " ")), r[1], r[2], "%.4f" % r[3], "%.4f" % r[4], "%+.4f" % r[5], "%d" % r[6], r[7], r[8])
    for (xx, ha), t in zip(cols, vals):
        ax.text(xx, y, t, ha=ha, va="center", fontsize=8.5)
    y -= 1
if cmp_rows:
    ax.plot([0, 1], [y + 0.5, y + 0.5], color="black", lw=0.6)
    y -= 0.25
    for (xx, ha), t in zip([(0.0, "left"), (0.15, "left"), (0.56, "right"), (1.0, "right")],
                           ("injection", "priors at", r"fast $-$ original, $-\ln L$", r"max $|$fast$-$original$|$/width")):
        ax.text(xx, y, t, ha=ha, va="center", fontsize=8.5)
    y -= 1
    for r in cmp_rows:
        for (xx, ha), t in zip([(0.0, "left"), (0.15, "left"), (0.56, "right"), (1.0, "right")],
                               (DATA_LABEL.get(r[0], r[0].replace("_", " ")), r[1], "%+.1e" % r[2], r[3])):
            ax.text(xx, y, t, ha=ha, va="center", fontsize=8.5)
        y -= 1
ax.plot([0, 1], [0, 0], color="black", lw=1.2)
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
fig.savefig(os.path.join(here, "plots", "injection_summary.pdf"))
fig.savefig(os.path.join(here, "plots", "injection_summary.png"))
plt.close(fig)

for r in rows:
    print("%-22s priors=%-7s %-8s  -lnL(truth) %.4f  best %.4f  best-truth %+.4f  evals %4d  %-12s  max|fit-truth| %s  succeeded=%s" % r)
for r in cmp_rows:
    print("%-22s priors=%-7s fast-original: d(-lnL) %+.2e, max |dparam| %s" % r)
