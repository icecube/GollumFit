#!/usr/bin/env python3
"""Compare two bench_llh.py outputs. Exit code 1 if any tolerance is violated.

Usage: python bench_compare.py reference.npz candidate.npz [--rtol-exp 1e-5] [--atol-llh 1e-3] [--atol-fit 1e-3]
"""
import argparse, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("reference"); ap.add_argument("candidate")
ap.add_argument("--rtol-exp", type=float, default=1e-5)
ap.add_argument("--atol-llh", type=float, default=5e-3)  # photospline float rounding gives ~1e-3 on LLH~2e4; 1 sigma is 0.5
ap.add_argument("--atol-fit", type=float, default=1e-3)
a = ap.parse_args()
R, C = np.load(a.reference), np.load(a.candidate)
ok = True

e_r, e_c = R["expectations"], C["expectations"]
mask = e_r > 0
rel = np.zeros_like(e_r); rel[mask] = np.abs(e_c[mask] - e_r[mask]) / e_r[mask]
absd = np.abs(e_c - e_r)
worst_rel = rel.max(); worst_abs = absd.max()
flag = "OK" if worst_rel <= a.rtol_exp else "FAIL"; ok &= flag == "OK"
print("expectations: max rel diff %.3e (max abs %.3e, totals %s)  [%s]" % (worst_rel, worst_abs,
      np.array2string(e_c.sum(axis=(1, 2, 3)) - e_r.sum(axis=(1, 2, 3)), precision=5), flag))

for key in ("llh_prior", "llh_noprior"):
    d = np.abs(C[key] - R[key]); flag = "OK" if d.max() <= a.atol_llh else "FAIL"; ok &= flag == "OK"
    print("%-12s max abs diff %.3e  (ref %s)  [%s]" % (key, d.max(), np.array2string(R[key], precision=4), flag))

if np.isfinite(R["fit_llh"]) and np.isfinite(C["fit_llh"]):
    d = np.abs(C["fit_params"] - R["fit_params"]); dl = abs(float(C["fit_llh"]) - float(R["fit_llh"]))
    flag = "OK" if (d.max() <= a.atol_fit and dl <= a.atol_llh) else "FAIL"; ok &= flag == "OK"
    print("fit: max |dparam| %.3e, |dLLH| %.3e, nEval %d -> %d  [%s]" % (d.max(), dl, R["fit_neval"], C["fit_neval"], flag))

def t(x, k):
    return float(np.min(x[k])) if np.size(x[k]) else float("nan")
print("timing: EvalLLH %.3fs -> %.3fs (x%.2f) | GetExpectation %.3fs -> %.3fs | per-gradient-eval %.3fs -> %.3fs (x%.2f) | threads %d -> %d" % (
    t(R, "t_eval"), t(C, "t_eval"), t(R, "t_eval") / t(C, "t_eval"),
    t(R, "t_exp"), t(C, "t_exp"),
    float(R["fit_wall"]) / max(1, int(R["fit_neval"])), float(C["fit_wall"]) / max(1, int(C["fit_neval"])),
    (float(R["fit_wall"]) / max(1, int(R["fit_neval"]))) / (float(C["fit_wall"]) / max(1, int(C["fit_neval"]))),
    int(R["threads"]), int(C["threads"])))
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
