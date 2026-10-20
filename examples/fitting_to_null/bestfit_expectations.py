#!/usr/bin/env python3
"""Expected bin contents and likelihood at the best-fit points of finished fits, for one build.

Run once per build, in that build's environment (the two builds cannot be imported together):
  source GollumFit_test/setup_head.sh
  python bestfit_expectations.py --out bestfit_head.npz head8=<headfit log>:0 head1=<headfit log>:1
  source GollumFit_test/setup_perf.sh
  python bestfit_expectations.py --out bestfit_fast.npz fast1=fitcmp_<job>.out fast8=fitcmp_<job>.out
Each positional argument is name=path[:index]; index selects the n-th fit in the log (default 0).
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("fits", nargs="+")
args = ap.parse_args()

import GollumFitPy as gf
g = bc.make_gollumfit(gf, threads=1)
g.SetData(bc.load_null_data())
flag, bound, priors = bc.make_fit_config(gf)
g.SetFitParametersFlag(flag); g.SetFitParametersBound(bound); g.SetFitParametersPriors(priors)
g.SetFitParametersSeed([bc.nominal_params(gf)])
g.ConstructLikelihoodProblem()

out = {"data": np.asarray(g.GetDataDistribution()), "module": gf.__file__}
for spec in args.fits:
    name, rest = spec.split("=", 1)
    path, _, idx = rest.partition(":")
    fit = bc.parse_fits(path)[int(idx or 0)]
    fp = bc.nominal_params(gf, fit["params"])
    out[name + "_expectation"] = np.asarray(g.GetExpectation(fp))
    out[name + "_llh_recomputed"] = g.EvalLLH(fp, True)
    out[name + "_params"] = bc.params_to_array(fp)
    out[name + "_llh"] = fit["llh"]
    out[name + "_neval"] = fit["neval"]
    out[name + "_wall"] = fit["wall"]
    print("[bestfit] %s: LLH %.6f (recomputed %.6f), nEval %d, wall %.1f s, total expectation %.3f"
          % (name, fit["llh"], out[name + "_llh_recomputed"], fit["neval"], fit["wall"], out[name + "_expectation"].sum()),
          file=sys.stderr, flush=True)
np.savez(args.out, **out)
print("[bestfit] wrote", args.out, file=sys.stderr)
