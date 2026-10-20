#!/usr/bin/env python3
"""Fit injected fake data (from generate_injection.py) with all 38 parameters free.

Same configuration and starting point as generate_likelihood.py (seed 100 throw from the priors),
but with the fake data and, optionally, the Gaussian prior centres moved to the injected values.

Usage: python fit_injection.py --data injection_data.npz --priors {truth,nominal}
       [--threads 1] [--change-tol 1e-20] [--grad-tol 1e-20]
Prints in the format read by bench_common.parse_fits, plus the likelihood at the injected point.
"""
import argparse, os, sys, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="injection_data.npz")
ap.add_argument("--priors", choices=("truth", "nominal"), required=True)
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--change-tol", type=float, default=1e-20)
ap.add_argument("--grad-tol", type=float, default=1e-20)
args = ap.parse_args()

d = np.load(args.data)
truth = dict(zip([str(n) for n in d["names"]], d["truth"]))

import GollumFitPy as gf
print("module", gf.__file__, flush=True)
g = bc.make_gollumfit(gf, threads=args.threads, change_tol=args.change_tol, grad_tol=args.grad_tol)
g.SetData(d["realization"])
flag, bound, priors = bc.make_fit_config(gf, free=bc.NAMES)
if args.priors == "truth":
    for name, (ptype, *_) in bc.SYST.items():
        if ptype == "Gaussian":
            setattr(priors, name + "Center", float(truth[name]))

# starting point: identical to generate_likelihood.py (np.random.seed(100), throw in this order)
np.random.seed(100)
seed = gf.FitParameters()
for name, (ptype, c, w, lo, hi) in bc.SYST.items():
    if ptype == "Gaussian":
        v = stats.truncnorm((lo - c) / w, (hi - c) / w, c, w).rvs(1)[0]
    else:
        v = np.random.uniform(lo, hi)
    setattr(seed, name, float(v))

g.SetFitParametersFlag(flag); g.SetFitParametersBound(bound); g.SetFitParametersPriors(priors)
g.SetFitParametersSeed([seed])
g.ConstructLikelihoodProblem()

llh_truth = g.EvalLLH(bc.nominal_params(gf, truth), True)
llh_start = g.EvalLLH(seed, True)
print("priors centred at:", args.priors, flush=True)
print("LLH at truth: %.9f" % llh_truth, flush=True)
print("LLH at start: %.9f" % llh_start, flush=True)
t0 = time.time()
r = g.MinLLH()
print("MinLLH wall time: %.1f s" % (time.time() - t0))
for name in bc.NAMES:
    print(name, getattr(r.params, name))
print("llh:", r.likelihood)
print("nEval:", r.nEval)
print("succeeded:", r.succeeded, flush=True)
