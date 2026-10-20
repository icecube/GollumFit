#!/usr/bin/env python3
"""Resident memory of a GollumFit likelihood session, stage by stage.

Usage: python mem_probe.py --compact FILE --threads N --label TEXT
Prints one JSON line: for each stage the current resident set (VmRSS) and the peak so far (VmHWM),
in MB, from /proc/self/status. Stages: python + numpy; module import; GollumFit constructor (compact
file, splines, simulation histogram); data + likelihood problem; first value evaluation (the fast
build builds its per-event cache here); four more value evaluations; a short fit (value + gradient
evaluations, 1 free parameter, loose tolerances).
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc


def mem():
    out = {}
    for line in open("/proc/self/status"):
        if line.startswith(("VmRSS:", "VmHWM:")):
            out[line.split(":")[0]] = int(line.split()[1]) / 1024.0
    return out["VmRSS"], out["VmHWM"]


ap = argparse.ArgumentParser()
ap.add_argument("--compact", required=True)
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--label", default="")
args = ap.parse_args()

stages = []
def record(name, **extra):
    rss, hwm = mem()
    stages.append(dict(stage=name, rss_mb=round(rss, 1), peak_mb=round(hwm, 1), **extra))

record("python + numpy")
import GollumFitPy as gf
record("module import")
t0 = time.time()
g = bc.make_gollumfit(gf, threads=args.threads, change_tol=1e-1, grad_tol=1e-1, compact=args.compact)
record("constructor (compact file, splines, histogram)", seconds=round(time.time() - t0, 1))
g.SetData(bc.load_null_data())
flag, bound, priors = bc.make_fit_config(gf, free=["convNorm"])
g.SetFitParametersFlag(flag); g.SetFitParametersBound(bound); g.SetFitParametersPriors(priors)
g.SetFitParametersSeed([bc.nominal_params(gf, {"convNorm": 1.15})])
g.ConstructLikelihoodProblem()
record("data + likelihood problem")
fp = bc.nominal_params(gf)
t0 = time.time(); g.EvalLLH(fp, True)
record("first value evaluation", seconds=round(time.time() - t0, 2))
t0 = time.time()
for _ in range(4):
    g.EvalLLH(fp, True)
record("4 more value evaluations", seconds=round((time.time() - t0) / 4, 3))
t0 = time.time(); r = g.MinLLH()
record("short fit (value + gradient)", seconds=round(time.time() - t0, 1), n_eval=int(r.nEval))

print("MEMJSON " + json.dumps(dict(label=args.label, module=gf.__file__, threads=args.threads,
                                    compact=os.path.basename(args.compact), stages=stages)), flush=True)
