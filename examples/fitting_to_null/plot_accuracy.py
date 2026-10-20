#!/usr/bin/env python3
"""Timing and accuracy figures: original build (committed HEAD) vs fast build, 1 and 8 threads.

Inputs (examples/fitting_to_null):
  bestfit_head.npz, bestfit_fast.npz      from bestfit_expectations.py (best fits of the four runs)
  bench_head.npz, bench_final.npz         from bench_llh.py --no-fit (both builds at five fixed parameter points)
  timing_value_<head|fast>_t<1|8>.npz     from bench_value_only.slurm (value-only EvalLLH timings)
Outputs (plots/, .pdf and .png):
  timing                    fit wall time, time per fit evaluation, time per value-only evaluation
  accuracy_llh              table of the four best-fit likelihood values
  accuracy_bestfit          best-fit parameters (pulls) and differences to the original 1-thread fit
  accuracy_expectations     best-fit track spectra and ratios to the original 1-thread fit
  accuracy_bin_differences  per-bin best-fit expectation, fast vs original at 1 and at 8 threads
  accuracy_same_parameters  per-bin expectation, fast vs original at identical parameters (code only)
Run with PYTHONNOUSERSITE=1 on fasrc (a broken PIL in ~/.local otherwise breaks matplotlib).
"""
import os, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter, MaxNLocator

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import bench_common as bc

plt.style.use(os.path.join(here, "figures.mplstyle"))
if subprocess.call(["which", "latex"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
    plt.rcParams.update({"text.usetex": False, "font.family": "serif", "mathtext.fontset": "cm"})
plt.rcParams.update({"axes.labelsize": 13, "xtick.labelsize": 11, "ytick.labelsize": 11})
os.makedirs(os.path.join(here, "plots"), exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(here, "plots", name + ".pdf"))
    fig.savefig(os.path.join(here, "plots", name + ".png"))
    plt.close(fig)


def pow10_scale(x):
    """Exponent e such that x / 10^e is of order 1 (used in axis labels instead of offset text)."""
    return int(np.floor(np.log10(x))) if x > 0 else 0


def fmt_duration(t):
    """Human-readable duration: '57 min 41 s', '1 h 3 min', '15.5 s', '124 ms'."""
    if t >= 3600:
        h, m = int(t // 3600), int(round((t % 3600) / 60))
        if m == 60:
            h, m = h + 1, 0
        return "%d h %d min" % (h, m)
    if t >= 60:
        m, s = int(t // 60), int(round(t % 60))
        if s == 60:
            m, s = m + 1, 0
        return "%d min %d s" % (m, s)
    if t >= 1:
        return "%.1f s" % t if t < 10 else "%.0f s" % t
    ms = 1e3 * t
    return "%.1f ms" % ms if ms < 10 else "%.0f ms" % ms


C_OLD, C_NEW = "#4C72B0", "#DD8452"
CASES = [  # key, label, color, linestyle, marker
    ("head1", "original, 1 thread", C_OLD, "-", "o"),
    ("head8", "original, 8 threads", C_OLD, "--", "^"),
    ("fast1", "fast, 1 thread", C_NEW, "-", "s"),
    ("fast8", "fast, 8 threads", C_NEW, "--", "D"),
]
LABEL = {k: l for k, l, *_ in CASES}
REF = "head1"
F = {}
for f in ("bestfit_head.npz", "bestfit_fast.npz"):
    d = np.load(os.path.join(here, f), allow_pickle=True)
    for k in d.files:
        F[k] = d[k]
cases = [c for c in CASES if c[0] + "_wall" in F]
missing = [c[0] for c in CASES if c[0] + "_wall" not in F]
if missing:
    print("warning: missing cases", missing)

# ======================= timing =======================
value_only = {}
for key, *_ in CASES:
    f = os.path.join(here, "timing_value_%s_t%s.npz" % (key[:-1], key[-1]))
    if os.path.exists(f):
        value_only[key] = float(np.mean(np.load(f)["t_eval"]))
quantities = ["fit", "eval"] + (["value"] if value_only else [])
fig, axes = plt.subplots(1, len(quantities), figsize=(4 * len(quantities), 3.2))
pos = {"head1": 0.0, "fast1": 1.0, "head8": 2.5, "fast8": 3.5}
DUR_TICKS = [(1e-3, "1 ms"), (1e-2, "10 ms"), (1e-1, "100 ms"), (1, "1 s"), (10, "10 s"),
             (60, "1 min"), (600, "10 min"), (3600, "1 h"), (36000, "10 h")]
for ax, quantity in zip(axes, quantities):
    vals = {}
    for key, label, color, ls, _ in cases:
        if quantity == "value" and key not in value_only:
            continue
        wall, neval = float(F[key + "_wall"]), int(F[key + "_neval"])
        v = wall if quantity == "fit" else (wall / neval if quantity == "eval" else value_only[key])
        vals[key] = v
        ax.bar(pos[key], v, width=0.8, color=color, alpha=0.85 if key.endswith("1") else 0.55,
               hatch=None if key.endswith("1") else "//", edgecolor=color, linewidth=0.8)
        ax.text(pos[key], v * 1.25, fmt_duration(v), ha="center", fontsize=8)
    for t in ("1", "8"):
        if "head" + t in vals and "fast" + t in vals:
            r = vals["head" + t] / vals["fast" + t]
            ax.text((pos["head" + t] + pos["fast" + t]) / 2, vals["head" + t] * 4.0, r"$\times%.0f$" % r,
                    ha="center", fontsize=10, color="0.25")
    ax.set_yscale("log")
    lo, hi = min(vals.values()) / 3, max(vals.values()) * 15
    ax.set_ylim(lo, hi)
    ticks = [(v, l) for v, l in DUR_TICKS if lo <= v <= hi]
    ax.yaxis.set_major_locator(FixedLocator([v for v, _ in ticks]))
    ax.yaxis.set_major_formatter(FixedFormatter([l for _, l in ticks]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="y", labelsize=9)
    ax.set_xticks([pos[k] for k, *_ in cases])
    ax.set_xticklabels([l.replace(", ", "\n") for _, l, *_ in cases], fontsize=8)
    ax.tick_params(axis="x", which="both", length=0)
    ax.set_ylabel({"fit": "38-parameter fit wall time",
                   "eval": "Time per fit evaluation\n(value + 38-parameter gradient)",
                   "value": "Time per likelihood evaluation\n(value only)"}[quantity], fontsize=10)
fig.tight_layout()
save(fig, "timing")

# ======================= best-fit likelihood values (table) =======================
def sci_tex(x):
    """'0' or '+1.0\times10^{-3}' style for LaTeX."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    return r"$%+.1f\times10^{%d}$" % (x / 10.0 ** e, e)

llh_ref = float(F[REF + "_llh"])
rows = []
for key, label, *_ in cases:
    llh = float(F[key + "_llh"])
    d = llh - llh_ref
    rows.append((label, "%.6f" % llh, "0 (baseline)" if key == REF else sci_tex(d)))
HEAD_H = 1.5  # header row height in row units (two-line header)
n = len(rows) + HEAD_H
fig, ax = plt.subplots(figsize=(4.6, 0.3 * n + 0.2))
ax.axis("off")
ax.set_xlim(0, 1)
ax.set_ylim(0, n)
cols = [(0.02, "left"), (0.60, "right"), (0.98, "right")]
header = ("fit", r"$-\ln L$ at best fit", "difference to\noriginal, 1 thread")
for (x, ha), txt in zip(cols, header):
    ax.text(x, n - HEAD_H / 2, txt, ha=ha, va="center", fontsize=10, linespacing=1.1)
for i, row in enumerate(rows):
    for (x, ha), txt in zip(cols, row):
        ax.text(x, len(rows) - 0.5 - i, txt, ha=ha, va="center", fontsize=10)
ax.plot([0, 1], [n, n], color="black", lw=1.2, clip_on=False)
ax.plot([0, 1], [len(rows), len(rows)], color="black", lw=0.6, clip_on=False)
ax.plot([0, 1], [0, 0], color="black", lw=1.2, clip_on=False)
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
save(fig, "accuracy_llh")

# ======================= best-fit parameters =======================
names = bc.NAMES
center = np.array([bc.SYST[n][1] for n in names])
width = np.array([bc.SYST[n][2] for n in names])
p_ref = F[REF + "_params"]
x = np.arange(len(names))
offsets = {"head1": -0.27, "head8": -0.09, "fast1": 0.09, "fast8": 0.27}
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True, gridspec_kw={"height_ratios": [1.2, 1], "hspace": 0.08})
dmax = 0.0
for key, label, color, ls, marker in cases:
    p = F[key + "_params"]
    style = dict(color=color, ms=3.5, ls="none", mfc=color if key.endswith("1") else "none", mew=1.0)
    ax1.plot(x + offsets[key], (p - center) / width, marker, label=label, **style)
    delta = (p - p_ref) / width
    dmax = max(dmax, np.abs(delta).max())
    dl = float(F[key + "_llh"]) - llh_ref
    ax2.plot(x + offsets[key], delta, marker,
             label="%s (baseline)" % label if key == REF else r"%s: $\Delta(-\ln L)$ = %+.1e" % (label, dl), **style)
ax1.axhline(0, color="gray", lw=0.8)
ax1.set_ylabel("(fit $-$ prior center)\n/ prior width", fontsize=10)
ax1.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=4, frameon=False)
ax2.axhline(0, color="gray", lw=0.8)
lim = max(0.05, 1.3 * dmax)
ax2.set_ylim(-lim, lim)
ax2.set_ylabel("(fit $-$ original, 1 thread)\n/ prior width", fontsize=10)
ax2.legend(fontsize=7, loc="lower left", ncol=2, frameon=False)
ax2.set_xticks(x)
ax2.set_xticklabels(names, rotation=90, fontsize=7)
for ax in (ax1, ax2):
    ax.tick_params(axis="x", which="both", length=0)
save(fig, "accuracy_bestfit")

# ======================= expected track spectra =======================
TRACK = 1  # topology index: 0 = shower, 1 = track
edges = np.logspace(np.log10(300), np.log10(1e5), 25)
cedges = np.linspace(-1, 0, 21)
data = F["data"]
ref = F[REF + "_expectation"]
ratios, rmax = [], 0.0
fig, axes = plt.subplots(2, 2, figsize=(8, 4.8), sharex="col", gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.06, "wspace": 0.3})
for col, (axis_sum, xedges) in enumerate(((0, edges), (1, cedges))):
    dd = data[TRACK].sum(axis=axis_sum)
    axes[0, col].fill_between(xedges, np.r_[dd, dd[-1]], step="post", color="0.88", lw=0, zorder=0,
                              label="fake data" if col == 0 else None)
    for key, label, color, ls, _ in cases:
        is_head = key.startswith("head")
        y = F[key + "_expectation"][TRACK].sum(axis=axis_sum)
        # original build drawn wide underneath, fast build narrow on top, so overlapping curves stay visible
        axes[0, col].step(xedges, np.r_[y, y[-1]], where="post", color=color, linestyle=ls,
                          lw=3.0 if is_head else 1.1, alpha=0.6 if is_head else 1.0, zorder=1 if is_head else 2,
                          label=label if col == 0 else None)
        if key != REF:
            r = y / ref[TRACK].sum(axis=axis_sum) - 1
            rmax = max(rmax, np.abs(r).max())
            ratios.append((col, xedges, r, color, ls, key))
ex = pow10_scale(rmax)
for col, xedges, r, color, lstyle, key in ratios:
    rs = r / 10.0 ** ex
    axes[1, col].step(xedges, np.r_[rs, rs[-1]], where="post", color=color, linestyle=lstyle,
                      lw=2.4 if key.startswith("head") else 1.2)
axes[0, 0].set_xscale("log")
for ax in axes[0]:
    ax.set_yscale("log")
axes[0, 0].set_ylabel("Expected track events")
axes[1, 0].set_xlabel("Reconstructed energy [GeV]")
axes[1, 1].set_xlabel(r"$\cos\theta_{\mathrm{reco}}$")
lim = 1.3 * rmax / 10.0 ** ex if rmax > 0 else 1.0
for ax in axes[1]:
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_ylim(-lim, lim)
    ax.set_ylabel(r"ratio $-$ 1 [$10^{%d}$]" % ex, fontsize=10)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
axes[0, 0].legend(fontsize=7.5, loc="lower left")
axes[0, 0].text(0.97, 0.95, "tracks\nlower panels: ratio to\noriginal, 1 thread",
                transform=axes[0, 0].transAxes, fontsize=7.5, va="top", ha="right")
for ax in axes.flat:
    ax.tick_params(which="minor", length=3)
save(fig, "accuracy_expectations")

# ======================= per-bin best-fit expectation: fast vs original =======================
pairs = [("head1", "fast1", "1 thread"), ("head8", "fast8", "8 threads")]
pairs = [p for p in pairs if p[0] + "_expectation" in F and p[1] + "_expectation" in F]
diffs = {}
for a, b, _ in pairs:
    ea, eb = F[a + "_expectation"], F[b + "_expectation"]
    m = ea > 0
    diffs[(a, b)] = (eb[m] - ea[m]) / ea[m]
span = max(np.abs(v).max() for v in diffs.values())
e2 = pow10_scale(span)
bins = np.linspace(-1.05, 1.05, 50) * span / 10.0 ** e2
fig, axes = plt.subplots(1, len(pairs), figsize=(4 * len(pairs), 3.2), sharey=True, gridspec_kw={"wspace": 0.12})
axes = np.atleast_1d(axes)
for ax, (a, b, tlabel) in zip(axes, pairs):
    d = diffs[(a, b)]
    ax.hist(d / 10.0 ** e2, bins=bins, histtype="step", lw=1.5, color=C_NEW,
            linestyle="-" if tlabel.startswith("1") else "--")
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_xlabel(r"fast/original $-$ 1 per bin [$10^{%d}$]" % e2, fontsize=11)
    ax.text(0.97, 0.95, "best fits, %s\n%d bins\nmedian $|\\Delta|$ %.0e\nmax $|\\Delta|$ %.1e"
            % (tlabel, d.size, np.median(np.abs(d)), np.abs(d).max()),
            transform=ax.transAxes, fontsize=8, va="top", ha="right")
    ax.xaxis.set_major_locator(FixedLocator([v for v in (-1.0, -0.5, 0.0, 0.5, 1.0) if abs(v) < bins[-1]]))
    ax.tick_params(which="minor", length=3)
axes[0].set_ylabel("Bins", fontsize=11)
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.2, top=0.96)
save(fig, "accuracy_bin_differences")

# ======================= per-bin expectation at identical parameters (code only) =======================
H = np.load(os.path.join(here, "bench_head.npz"))
N = np.load(os.path.join(here, "bench_final.npz"))
rel = []
for k in range(len(bc.TEST_POINTS)):
    a_, b_ = H["expectations"][k], N["expectations"][k]
    m = a_ > 0
    rel.append((b_[m] - a_[m]) / a_[m])
rel = np.concatenate(rel)
e1 = pow10_scale(np.abs(rel).max())
fig, ax = plt.subplots(figsize=(5, 3.2))
ax.hist(rel / 10.0 ** e1, bins=50, histtype="step", lw=1.5, color=C_NEW)
ax.axvline(0, color="gray", lw=0.8)
ax.set_xlabel(r"fast/original $-$ 1 per bin [$10^{%d}$]" % e1, fontsize=11)
ax.set_ylabel("Bins", fontsize=11)
ax.text(0.97, 0.95, "identical parameters\n(5 points, %d bins)\nmax $|\\Delta|$ %.1e" % (rel.size, np.abs(rel).max()),
        transform=ax.transAxes, fontsize=8, ha="right", va="top")
ax.tick_params(which="minor", length=3)
fig.tight_layout()
save(fig, "accuracy_same_parameters")

print("wrote", sorted(os.listdir(os.path.join(here, "plots"))))
