#!/usr/bin/env python3
"""Fake data at randomly perturbed nuisance parameters, for an injection-recovery test.

Draws all 38 nuisance parameters from the priors used in the fit, truncated to the fit bounds
(uniform for astroPivot), with the configuration of generate_likelihood.py, and writes the expected
events at that point (Asimov data, no Poisson fluctuation) in the format of nullexpectation.npz.
By default the 16 hadronic/cosmic-ray and the 9 ice-gradient parameters are drawn jointly from
their correlated Gaussian priors (the ones in resources/correlation_matrices), as the fit uses them;
--independent draws every parameter separately and ignores the correlations (the flux correlation
matrix is nearly degenerate, so such a point can be very unlikely under the fit's prior).
Run it with the original build (setup_head.sh) so the fake data come from the reference code.

Usage: python generate_injection.py [--seed 2027] [--independent] [--out injection_data.npz]
Output keys: realization (N x 4), truth (38, the float32-rounded injected values), names, seed.
"""
import argparse, os, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc

ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int, default=2027)
ap.add_argument("--independent", action="store_true", help="ignore the prior correlations")
ap.add_argument("--out", default="injection_data.npz")
args = ap.parse_args()

rs = np.random.RandomState(args.seed)
truth = {}
for name, (ptype, c, w, lo, hi) in bc.SYST.items():
    if ptype == "Gaussian":
        truth[name] = float(stats.truncnorm((lo - c) / w, (hi - c) / w, c, w).rvs(random_state=rs))
    else:
        truth[name] = float(rs.uniform(lo, hi))

FLUX = ["hadronicHEkp", "hadronicHEkm", "hadronicVHE1pip", "hadronicVHE1pim", "hadronicVHE3kp", "hadronicVHE3km",
        "hadronicVHE3pip", "hadronicVHE3pim", "hadronicVHE3p", "hadronicVHE3n",
        "cosmicRay1", "cosmicRay2", "cosmicRay3", "cosmicRay4", "cosmicRay5", "cosmicRay6"]
ICE = ["icegrad%d" % i for i in range(9)]
root = bc.repo_root()
blocks = [(FLUX, "resources/correlation_matrices/flux_correlations_new_ddmnodeis.npy"),
          (ICE, "resources/correlation_matrices/icegrad_correlations.npy")]
if not args.independent:
    for block, f in blocks:
        corr = np.load(os.path.join(root, f))
        c = np.array([bc.SYST[n][1] for n in block]); w = np.array([bc.SYST[n][2] for n in block])
        lo = np.array([bc.SYST[n][3] for n in block]); hi = np.array([bc.SYST[n][4] for n in block])
        cov = corr * np.outer(w, w)
        for attempt in range(100000):  # truncate to the fit bounds by rejection
            x = rs.multivariate_normal(c, cov)
            if np.all(x > lo) and np.all(x < hi):
                break
        for n, v in zip(block, x):
            truth[n] = float(v)
for block, f in blocks:
    corr = np.load(os.path.join(root, f))
    z = np.array([(truth[n] - bc.SYST[n][1]) / bc.SYST[n][2] for n in block])
    print("[inject] %s block: prior chi2 %.1f for %d parameters (correlated prior)" % ("flux" if block is FLUX else "ice", z @ np.linalg.solve(corr, z), len(z)), file=sys.stderr)

import GollumFitPy as gf
g = bc.make_gollumfit(gf, threads=1)
fp = bc.nominal_params(gf, truth)
truth_arr = bc.params_to_array(fp)  # what the code actually sees (float32)
events = np.asarray(g.GetExpectationEvents(fp))
np.savez(args.out, realization=events, truth=truth_arr, names=np.array(bc.NAMES), seed=args.seed)
print("[inject] module %s" % gf.__file__, file=sys.stderr)
print("[inject] seed %d, %d events, total expectation %.3f -> %s" % (args.seed, len(events), events[:, 3].sum(), args.out), file=sys.stderr)
for n, v in zip(bc.NAMES, truth_arr):
    c, w = bc.SYST[n][1], bc.SYST[n][2]
    print("[inject] %-26s %10.5f  (%+.2f prior widths from nominal)" % (n, v, (v - c) / w), file=sys.stderr)
