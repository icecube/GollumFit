#!/usr/bin/env python3
"""Memory of parallel likelihood workers, as an MCMC pool would use them.

Mode 'own':  every worker process builds its own GollumFit object (multiprocessing initializer).
Mode 'fork': the parent builds one object, evaluates it once (which builds the lazily created
             per-event table and bin cache), then forks the workers; they share its memory
             copy-on-write and only read it.
Each worker evaluates the likelihood at random points around nominal, then reports its memory
from /proc/self/smaps_rollup. PSS splits shared pages between the processes that share them, so
the sum of PSS over all processes is the physical memory the pool really uses.
Usage: python mcmc_fork_memory.py --mode fork|own --workers 16
"""
import argparse, multiprocessing as mp, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc
import GollumFitPy as gf

G = None


def smaps():
    out = {}
    for line in open("/proc/self/smaps_rollup"):
        parts = line.split()
        if parts and parts[0].rstrip(":") in ("Rss", "Pss", "Private_Dirty", "Shared_Clean", "Shared_Dirty"):
            out[parts[0].rstrip(":")] = int(parts[1]) / 1024.0
    return out


def build():
    global G
    G = bc.make_gollumfit(gf, threads=1)  # forked workers must evaluate single-threaded (see README)
    G.SetData(bc.load_null_data())
    flag, bound, priors = bc.make_fit_config(gf)
    G.SetFitParametersFlag(flag); G.SetFitParametersBound(bound); G.SetFitParametersPriors(priors)
    G.SetFitParametersSeed([bc.nominal_params(gf)])
    G.ConstructLikelihoodProblem()


def evaluate(task):
    seed, n = task
    rng = np.random.default_rng(seed)
    nominal = np.array([bc.SYST[k][1] for k in bc.NAMES]); width = np.array([bc.SYST[k][2] for k in bc.NAMES])
    fp = gf.FitParameters()
    for _ in range(n):
        for name, v in zip(bc.NAMES, nominal + 0.02 * width * rng.standard_normal(len(nominal))):
            setattr(fp, name, float(v))
        G.EvalLLH(fp, True)
    check = G.EvalLLH(bc.nominal_params(gf), True)
    time.sleep(1.0)  # keep every worker busy at once so all of them are measured
    return os.getpid(), smaps(), check


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("fork", "own"), required=True)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    ctx = mp.get_context("fork")
    t0 = time.time()
    if args.mode == "fork":
        build()
        ref = G.EvalLLH(bc.nominal_params(gf), True)  # warm-up: builds the per-event table in the parent
        pool = ctx.Pool(args.workers)
    else:
        ref = None
        pool = ctx.Pool(args.workers, initializer=build)
    results = pool.map(evaluate, [(i, 20) for i in range(args.workers)], chunksize=1)
    parent = smaps()
    pool.close(); pool.join()
    per = {}
    for pid, m, check in results:
        per[pid] = m
        if ref is None:
            ref = check
        assert abs(check - ref) < 1e-9, "worker likelihood differs from the reference"
    rss = [m["Rss"] for m in per.values()]; pss = [m["Pss"] for m in per.values()]
    priv = [m["Private_Dirty"] for m in per.values()]
    print("[pool] mode %-4s workers %d (distinct %d) | parent RSS %.0f MB PSS %.0f MB | worker RSS mean %.0f MB, "
          "PSS mean %.0f MB, private-dirty mean %.0f MB | total PSS (parent + workers) %.0f MB | "
          "naive sum of RSS %.0f MB | LLH check %.6f | %.0f s"
          % (args.mode, args.workers, len(per), parent["Rss"], parent["Pss"], np.mean(rss), np.mean(pss),
             np.mean(priv), parent["Pss"] + sum(pss), parent["Rss"] + sum(rss), ref, time.time() - t0), flush=True)
