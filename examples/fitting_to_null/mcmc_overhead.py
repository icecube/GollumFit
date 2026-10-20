#!/usr/bin/env python3
"""Cost of a Python-driven random-walk Metropolis step vs a bare EvalLLH call (one thread).

Usage: python mcmc_overhead.py N_STEPS
Measures: (1) setting all 38 FitParameters attributes from a numpy vector, (2) EvalLLH alone,
(3) full Metropolis steps (propose, set parameters, EvalLLH with priors, accept/reject).
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc
import GollumFitPy as gf

n_steps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
g = bc.make_gollumfit(gf, threads=1)
g.SetData(bc.load_null_data())
flag, bound, priors = bc.make_fit_config(gf)
g.SetFitParametersFlag(flag); g.SetFitParametersBound(bound); g.SetFitParametersPriors(priors)
g.SetFitParametersSeed([bc.nominal_params(gf)])
g.ConstructLikelihoodProblem()

names = bc.NAMES
x = np.array([bc.SYST[n][1] for n in names])
step = 0.01 * np.array([bc.SYST[n][2] for n in names])
fp = gf.FitParameters()

def set_params(v):
    for n, val in zip(names, v):
        setattr(fp, n, float(val))

# (1) parameter transfer only
t0 = time.perf_counter()
for _ in range(1000):
    set_params(x)
t_set = (time.perf_counter() - t0) / 1000

# (2) bare likelihood evaluation
set_params(x)
g.EvalLLH(fp, True)  # warm-up (builds caches)
reps = max(3, n_steps // 5)
t0 = time.perf_counter()
for _ in range(reps):
    g.EvalLLH(fp, True)
t_eval = (time.perf_counter() - t0) / reps

# (3) Metropolis steps
rng = np.random.default_rng(1)
set_params(x); lp = -g.EvalLLH(fp, True)
acc = 0
t0 = time.perf_counter()
for _ in range(n_steps):
    y = x + step * rng.standard_normal(len(x))
    set_params(y)
    lq = -g.EvalLLH(fp, True)
    if np.log(rng.random()) < lq - lp:
        x, lp, acc = y, lq, acc + 1
t_step = (time.perf_counter() - t0) / n_steps
print("[mcmc] %s: set 38 params %.1f us | EvalLLH %.4f s | Metropolis step %.4f s | Python overhead %.2f %% | acceptance %d/%d"
      % (gf.__file__.split("/")[6], 1e6 * t_set, t_eval, t_step, 100 * (t_step - t_eval) / t_eval, acc, n_steps), flush=True)
