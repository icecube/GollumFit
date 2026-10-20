#!/usr/bin/env python3
"""Cross-check the fast weighters against the legacy weighter chain at several parameter points.

Per-event differences of order 1e-5 (value) and 1e-4 (gradient) are expected: photospline evaluates the
3-D tables in float arithmetic, so the reference itself carries that rounding. The histogram-level
comparison in bench_compare.py is the strict gate."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_common as bc
import GollumFitPy as gf

g = bc.make_gollumfit(gf, threads=1)
worst = 0.0
for i, ov in enumerate(bc.TEST_POINTS):
    fp = bc.nominal_params(gf, ov)
    t0 = time.time()
    d = g.CheckFastWeighter(fp)
    print("[check] point %d: max rel diff %.3e (%.1fs)" % (i, d, time.time() - t0), flush=True)
    worst = max(worst, d)
print("[check] WORST %.3e -> %s" % (worst, "PASS" if worst < 2e-4 else "FAIL"), flush=True)
sys.exit(0 if worst < 2e-4 else 1)
