#!/usr/bin/env python3
"""Benchmark + regression snapshot for the GollumFit likelihood.

Times EvalLLH, GetExpectation and a short MinLLH, and stores expectation
histograms / LLH values at fixed parameter points so that bench_compare.py can
check that an optimised build reproduces the reference numerically.

Usage:  python bench_llh.py [--threads N] [--out FILE.npz] [--no-fit]
                            [--free convNorm,astroNorm] [--change-tol X --grad-tol Y]
"""
import argparse, os, sys, time, resource
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--out", default="bench_result.npz")
ap.add_argument("--no-fit", action="store_true")
ap.add_argument("--free", default="convNorm,astroNorm")
ap.add_argument("--change-tol", type=float, default=1e-20)
ap.add_argument("--grad-tol", type=float, default=1e-20)
ap.add_argument("--repeat", type=int, default=3)
args = ap.parse_args()

def log(*a):
    print("[bench]", *a, file=sys.stderr, flush=True)

t0 = time.time()
import GollumFitPy as gf
log("import GollumFitPy from", gf.__file__, "%.2fs" % (time.time() - t0))

t0 = time.time()
g = bc.make_gollumfit(gf, threads=args.threads, change_tol=args.change_tol, grad_tol=args.grad_tol)
t_load = time.time() - t0
data = bc.load_null_data()
g.SetData(data)
free = [s for s in args.free.split(",") if s]
flag, bound, priors = bc.make_fit_config(gf, free=free)
seed = bc.nominal_params(gf, {"convNorm": 1.15, "astroNorm": 1.0})
g.SetFitParametersFlag(flag); g.SetFitParametersBound(bound); g.SetFitParametersPriors(priors)
g.SetFitParametersSeed([seed])
g.ConstructLikelihoodProblem()
log("construct %.2fs, data rows %d, RSS %.0f MB" % (t_load, len(data), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

# --- numerical snapshot at fixed points ---
expectations, llh_prior, llh_noprior = [], [], []
for ov in bc.TEST_POINTS:
    fp = bc.nominal_params(gf, ov)
    expectations.append(np.asarray(g.GetExpectation(fp)))
    llh_prior.append(g.EvalLLH(fp, True))
    llh_noprior.append(g.EvalLLH(fp, False))
expectations = np.array(expectations)
log("snapshot: total expectation at nominal %.4f, LLH(nominal) %.6f" % (expectations[0].sum(), llh_prior[0]))

# --- timings ---
fp0 = bc.nominal_params(gf)
t_eval = []
for _ in range(args.repeat):
    t0 = time.time(); g.EvalLLH(fp0, True); t_eval.append(time.time() - t0)
t_exp = []
for _ in range(max(1, args.repeat - 1)):
    t0 = time.time(); g.GetExpectation(fp0); t_exp.append(time.time() - t0)
log("EvalLLH(double): min %.3fs (threads=%d)   GetExpectation: min %.3fs" % (min(t_eval), args.threads, min(t_exp)))

fit = dict(params=np.full(len(bc.NAMES), np.nan), llh=np.nan, neval=0, wall=np.nan, cpu=np.nan)
if not args.no_fit:
    c0, w0 = time.process_time(), time.time()
    r = g.MinLLH()
    fit["wall"] = time.time() - w0; fit["cpu"] = time.process_time() - c0
    fit["params"] = bc.params_to_array(r.params); fit["llh"] = r.likelihood; fit["neval"] = int(r.nEval)
    log("MinLLH free=%s: LLH %.6f nEval %d wall %.1fs cpu %.1fs -> %.3fs per value+gradient evaluation"
        % (",".join(free), r.likelihood, r.nEval, fit["wall"], fit["cpu"], fit["wall"] / max(1, r.nEval)))

np.savez(args.out,
         names=np.array(bc.NAMES), threads=args.threads,
         expectations=expectations, llh_prior=np.array(llh_prior), llh_noprior=np.array(llh_noprior),
         t_load=t_load, t_eval=np.array(t_eval), t_exp=np.array(t_exp),
         fit_params=fit["params"], fit_llh=fit["llh"], fit_neval=fit["neval"], fit_wall=fit["wall"], fit_cpu=fit["cpu"],
         free=np.array(free), change_tol=args.change_tol, grad_tol=args.grad_tol)
log("wrote", args.out)
