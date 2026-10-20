# Cookbook: porting the fast-likelihood speed-up to another GollumFit install

Audience: an agent asked to bring the GollumFit likelihood speed-up (about 25x per evaluation,
about 30x per fit) into a different GollumFit checkout, branch or analysis, and to prove that it
gives the same physics result. Read `README_PERFORMANCE.md` (next to this directory) first for
what the changes are and why they are exact; this file is the procedure.

Reference implementation: GollumFit branch `perf/fast-likelihood` (based on main 620eec7) with
PhysTools branch `perf/likelihood-engine`, commit `ced3adf` (based on 8ce44a7). On fasrc the working
copy is `/n/home01/awen/meows_home/GollumFit_test/GollumFit-perf`.

## 0. Ground rules

- Do not commit, push or open pull requests unless the user asks. Prepare changes; the user
  reviews and commits.
- Never modify or rebuild the user's existing install or install prefix. Work in a git worktree
  with its own install prefix, and build a separate untouched baseline the same way.
- Every stage has a numerical gate (section 5). Do not move on, and do not report a speed-up,
  until the gate passes. Report failures as failures.
- Always test with at least two evaluation threads as well as one (pitfall 2).

## 1. What is in this directory

| path | contents | applies cleanly to |
|---|---|---|
| `patches/01-phystools-likelihood-engine.patch` | PhysTools: persistent thread pool, bin-pair cache with load-balanced tasks, running Kahan sums, no per-evaluation print, `histogram::binRef`, `BFGS_Function` by reference | PhysTools 8ce44a7 (`vendor/PhysTools`) |
| `patches/02-gollumfit-fast-likelihood.patch` | GollumFit core: `ReducedSplines.h`, `SystematicsCache.h`, `FastWeighter.h` (new), `WeighterMaker` returning the fast weighters (old chain kept as `LegacyWeighter`), `std::vector<Event>` containers, `RegisterSimulationWithWeighter`, `CheckFastWeighter`, `SetData` fix, `FitResult` init, one-block compact read, build flags and header dependencies in `configure`, optional Diver, GIL release, test root path | GollumFit 620eec7 (main) |
| `patches/03-examples-benchmarks-docs.patch` | example script options (threads/tolerances from env, timing print), benchmark and plotting scripts, `figures.mplstyle`, `README_PERFORMANCE.md` | GollumFit 620eec7 |
| `slurm/build_and_fit_head.slurm`, `../examples/fitting_to_null/bench_*.slurm` | the baseline build and the benchmark jobs used on fasrc; paths and partitions are fasrc-specific | adapt by hand |

Both GollumFit patches were checked with `git apply --check` on a clean `git archive 620eec7`, and
the PhysTools patch on a clean `git archive 8ce44a7`. Patch 02 does **not** apply to
`inelasticity_dev` (section 4.4).

The patches are snapshots of the two branches. If the branches move on, regenerate them:
```
git -C <PhysTools repo> diff 8ce44a7 perf/likelihood-engine > patches/01-phystools-likelihood-engine.patch
git diff 620eec7 perf/fast-likelihood -- configure include src python tests > patches/02-gollumfit-fast-likelihood.patch
git diff 620eec7 perf/fast-likelihood -- examples docs README_PERFORMANCE.md > patches/03-examples-benchmarks-docs.patch
```

## 2. Assess the target install

Run these and write down the answers before changing anything.

```
cd <target GollumFit checkout>
git log --oneline -1; git status --short | head          # base commit; uncommitted changes?
git merge-base --is-ancestor 620eec7 HEAD && echo "descends from 620eec7"
git diff --stat 620eec7 HEAD -- include/analysisWeighting.h include/Event.h include/GollumFit.h \
    src/GollumFit.cpp include/compactIO.h src/compactIO.cpp include/GollumParameters.h \
    include/FastMode.h python/GollumFitPy.cpp python/setup.py configure
git ls-tree HEAD vendor/PhysTools                         # PhysTools commit (8ce44a7?)
for p in 02-gollumfit-fast-likelihood 03-examples-benchmarks-docs; do
  git apply --check <cookbook>/patches/$p.patch && echo "$p applies"; done
```

Then check the three things the fast weighter depends on:

1. **The weight formula.** Read `sterile::WeighterMaker::operator()` in
   `include/analysisWeighting.h` (in the reference it is renamed `LegacyWeighter`). The fast
   weighter reproduces exactly this structure:
   ```
   w = convNorm * (conv + prompt + astro) * icegrad          (enableTotalNorm; else convNorm*conv + prompt + astro)
   conv   = HoleIce_c DOMEff_c Att_c (1+adu*ad)(1+klu*kl) (convFlux + sum h_i had_i + sum cr_j cr_j)
   prompt = promptNorm HoleIce_p DOMEff_p Att_p promptFlux
   astro  = astroNorm HoleIce_a DOMEff_a Att_a astroFlux nn tilt
   icegrad = prod_k (1 + g_k icegrad_k)
   ```
   with the 38 parameters in the order of `GollumFit::ConvertFitParameters`. If the target's
   formula, parameter list or parameter order differs in any way, the closed-form code in
   `FastWeighterCore::evaluate` must be changed to match (section 4.3).
2. **The spline tables.** Every DOM-efficiency, hole-ice and attenuation table must be
   3-dimensional with the nuisance parameter on the last axis, order at most 7, and at most 256
   coefficients along the parameter axis:
   ```
   python -c "import photospline as ps, glob
   for f in sorted(glob.glob('<spline dir>/*.fits')):
       t = ps.SplineTable(f); print(f.split('/')[-1], t.ndim, t.order, [len(k) for k in t.knots], t.extents)"
   ```
   The event coordinates are fixed in `SystematicsCache::build`: DOM efficiency and hole ice
   use `(log10 e.energy, cos e.zenith)` (reconstructed), attenuation uses
   `(log10 e.primaryEnergy, cos e.primaryZenith)` (true), and attenuation is skipped for
   `cos e.primaryZenith > 0.1`, as in the original weighters. The nominal values are
   `kDOMEfficiencyNominal = 1.27` and `kHoleIceNominal = -1.0`, hard-coded like the original
   weighters' nominal evaluations. Check that the target uses the same.
3. **Topologies and histograms.** Tables are looked up by `(FluxComponent, e.topology)` (DOM
   efficiency, hole ice) and `(FluxComponent, e.primaryType)` (attenuation). Check how the
   target assigns `e.topology` and which keys its spline maps contain. The PhysTools patch is
   generic over the histogram tuple (any number of histograms of any dimension).

Decide the route: **A** if patches 02 and 03 apply (target is main or close to it); **B**
(manual port, section 4.4) otherwise.

## 3. Isolated workspace and baseline

Use two worktrees with their own install prefixes: one for the unmodified baseline at the
target commit, one for the port. Never reuse the user's prefix.

```
S=<directory holding the checkout and its install prefix>      # e.g. GollumFit_test
git -C $S/<checkout> worktree add --detach $S/<name>-base <target commit>
git -C $S/<checkout> worktree add -b perf/fast-likelihood $S/<name>-perf <target commit>
cp -a $S/<prefix> $S/<prefix>_base; cp -a $S/<prefix> $S/<prefix>_perf
```

For each copy:
- Rewrite the pkg-config files to the new prefix:
  `sed -i "s|$S/<prefix>\b|$S/<prefix>_perf|g" $S/<prefix>_perf/lib/pkgconfig/*.pc`.
- Delete stale GollumFitPy eggs in `lib/python3.*/site-packages` (any `UNKNOWN-*`,
  `GollumFitPy-0.0.0-*` and old-version directories). `easy-install.pth` may list them.
- Write the environment script by copying the user's and changing the `GOLLUMBUILDPATH` line
  to the new prefix **by hand**, then `grep -n GOLLUMBUILDPATH` it (pitfall 5).
- The submodules in a fresh worktree are empty. Copy the vendored PhysTools working tree from
  the user's checkout (`cp -a <checkout>/vendor/PhysTools/. <worktree>/vendor/PhysTools/`),
  remove its `.git` file, and in the port worktree turn it into a standalone repo at the same
  commit so its changes can be diffed:
  `git init -q; git fetch -q <checkout>/vendor/PhysTools HEAD; git update-ref refs/heads/main FETCH_HEAD; git symbolic-ref HEAD refs/heads/main; git reset -q`.
- In the port prefix, replace `include/PhysTools` by a symlink to
  `<port worktree>/vendor/PhysTools/PhysTools` so header edits are picked up.
- Copy the data the example needs (compact FastMC file, fake data) into both worktrees, and
  check that they belong together: the number of rows in the fake data equals the number of
  simulated events in the compact file (pitfall 11).

Build both the same way (configure, make, make install, `python setup.py install --prefix=$PREFIX`
in `python/`), and print `GollumFitPy.__file__` in every job so the build actually used is on
record (pitfall 7). If the Python build fails on a missing `diver` library, make Diver
optional as in patch 02 (`python/setup.py`, `python/GollumFitPy.cpp`); change nothing else in
the baseline.

Baseline records, made with the baseline build:
- `bench_llh.py --no-fit --out bench_base.npz`: bin expectations and likelihoods at five
  parameter points. `bench_common.py` holds the example configuration (spline directories,
  binning, priors, test points); edit it if the target's example differs.
- One full fit with the example script, on a compute node, at 1 and at 8 threads.

## 4. Apply the changes in stages

Rebuild with `make clean && make && make install` and reinstall the Python module after every
stage that touches a header (pitfall 1).

### 4.1 PhysTools engine (patch 01)

```
cd <port worktree>/vendor/PhysTools && git apply <cookbook>/patches/01-phystools-likelihood-engine.patch
```
Gate: `bench_compare.py bench_base.npz <new>.npz` passes with differences at the 1e-10 level
(the engine changes only the summation order), and a two-parameter fit (`bench_llh.py`) gives
the same optimum. Run it with `--threads 1` and `--threads 4`.

### 4.2 Build flags and header dependencies

In `configure` (which writes `Makefile`): `CXXFLAGS += -std=c++17 -O3 -g -DNDEBUG -fno-math-errno`,
and add the PhysTools headers and the three new headers to the `GollumFit.o` dependency line.
`-march=haswell` (AVX2/FMA) was used in the reference build. It is not portable: it fails on
non-x86 machines and pre-2013 CPUs. Add it only if the user's machines support it, and say
so. GCC 9 has no `-march=x86-64-v3`. Do not use `-ffast-math` (the code relies on NaN checks
and sentinel values). Mirror the flags in `python/setup.py` `extra_compile_args`.

### 4.3 Reduced splines and the fast weighter (rest of patch 02)

Route A: `git apply <cookbook>/patches/02-gollumfit-fast-likelihood.patch` (and 03).

What the patch does, so it can be adapted:
1. New headers `include/ReducedSplines.h`, `include/SystematicsCache.h`, `include/FastWeighter.h`
   (self-contained; add them to the installed header list in `configure`).
2. `include/analysisWeighting.h`: include the new headers; give `sterile::WeighterMaker` a
   `shared_ptr<const SystematicsCache> cache_` and an `eventRanges_` list (copied in its copy
   and move operations; reset by the three `Set*Splines` methods); add
   `RegisterEventRanges(...)`, `GetCache()`; rename the templated `operator()` to
   `LegacyWeighter<DataType>(params)` and add two non-template overloads of `operator()` for
   `std::vector<double>` (returns `FastWeighter`) and `std::vector<FD<38>>` (returns
   `FastGradWeighter`).
3. Event storage: `std::deque<Event>` to `std::vector<Event>` for the simulation, meta-event and
   data containers (`GollumFit.h`, `FastMode.h`, `compactIO.h/.cpp`). The cache finds events by
   address, so a container must not grow after it is registered (pitfall 20).
4. `GollumFit::RegisterSimulationWithWeighter()`: registers `mainSimulation_` and
   `metaEvents_` with the weighter maker. Call it at the end of `ConstructSimulationHistogram`,
   at the end of `ConstructFastMode`, and in `SetSteeringParams` after `DFWM` is re-created.
5. `GollumFit::CheckFastWeighter(FitParameters)` plus its Python binding: compares every event's
   weight and 38 derivatives between `FastWeighter`/`FastGradWeighter` and `LegacyWeighter`.
6. Fixes: `SetData` resets `prob_` and `likelihood_problem_constructed_` (the problem holds
   references to the old data events); `FitResult` members initialised; compact file read in
   one block; GIL released (`py::call_guard<py::gil_scoped_release>()`) on the heavy bindings.
7. Defaults: `SteeringParams::change_tol` and `grad_tol` 1e-20 to 1e-6. This is a behaviour
   change. Ask the user before keeping it, and validate it on the final build (pitfall 18).

**If the weight formula differs from the reference**, edit `FastWeighterCore::evaluate`:
- Compute every factor once, then `w`, then write `g[i] = dw/dp_i` for all 38 parameters. A
  product `w = prod f_k` has `dw/dp = w * (d f_k/dp) / f_k`; use the explicit partial products
  (as for the ice gradients' prefix/suffix products) instead of dividing, so that factors equal
  to zero are safe.
- Spline factors: the value is `ReducedAxisSpline::eval(coef, basis)` and its derivative is
  `evalDeriv`. A DOM-efficiency or hole-ice factor `10^(r(p) - r0)` has derivative
  `ln10 * factor * r'(p)`.
- Keep `LegacyWeighter` as the definition of the weight and update it in the same change;
  `CheckFastWeighter` is the proof that the two agree, value and gradient.
- If a new spline needs other event coordinates or another nominal value, change
  `SystematicsCache::build` and its constants; if a table is not 3-D or the parameter is not
  its last axis, generalise `ReducedAxisSpline` (contract all axes except the kept one). Never
  silently fall back to approximate behaviour.

### 4.4 Route B: manual port to a diverged branch (e.g. `inelasticity_dev`)

Checked on `origin/inelasticity_dev`: it has the same PhysTools commit (patch 01 applies) and
an `analysisWeighting.h` that differs from main by one include, but `Event.h` (new
`inelasticity` field, new `readFile` signature), `GollumParameters.h`, `GollumFit.h`,
`src/GollumFit.cpp` (two histogram types: starting with an inelasticity axis, throughgoing;
two likelihood types), `compactIO` (metadata) and the Python bindings differ, so patch 02
fails on `configure`, `GollumParameters.h`, `compactIO.h` and `GollumFitPy.cpp`.

Procedure:
1. Apply patch 01 to the vendored PhysTools; the bin-pair cache is built per histogram in the
   tuple, so two histogram types are fine.
2. Copy the three new headers unchanged.
3. Apply the `analysisWeighting.h` part of patch 02 by hand (`git apply --include=include/analysisWeighting.h`
   usually works as the file barely differs).
4. Re-do items 3 to 6 of section 4.3 by hand in the branch's own versions of the files, keeping
   the branch's new logic (for example its compact-file metadata and readers).
5. Check how starting and throughgoing events set `e.topology` and that the loaded spline maps
   have keys for every topology that occurs; the cache throws if a corrected component has no
   table for an event's topology.
6. Adapt `bench_common.py` to the branch's steering parameters and data, then run all gates.

## 5. Gates (all must pass)

| check | command | pass if |
|---|---|---|
| every event, value and gradient | `python check_fast_weighter.py` | max relative difference below 2e-4 at all test points (typical: 1e-5 value, 7e-5 gradient: photospline's single-precision rounding) |
| bins and likelihood at fixed parameters | `bench_llh.py --no-fit`, then `bench_compare.py bench_base.npz <new>.npz` | expectations within 1e-5 relative, likelihood within 5e-3 (reference: 2.5e-6 and 1.4e-3) |
| unit tests | `pytest tests/test_gollumfit.py` | all pass (reference: 42/42) |
| threading | `bench_llh.py --threads 1`, `--threads 2`, `--threads 4` | identical numbers, no crash |
| full fit | example fit with both builds, same node type, 1 and 8 threads | likelihood within about 1e-3 (one sigma is 0.5); parameters within a few hundredths of a prior width except along flat directions |
| determinism | repeat one fit | same likelihood to every digit and same number of evaluations |

## 6. Benchmark protocol

- Compare against the baseline build made in section 3 from the clean target commit, with its
  original flags, never against an install of unknown provenance (pitfalls 10 and 12).
- Run all timing jobs on compute nodes of one CPU type. On fasrc the `test` partition starts
  immediately (12 h limit) when the lab partition is full; record `hostname` and `lscpu` model.
- Use the same script, seed, data and tolerances for both builds. Use the patched
  `generate_likelihood.py` (env `GOLLUM_THREADS`, `GOLLUM_CHANGE_TOL`, `GOLLUM_GRAD_TOL`; prints
  `MinLLH wall time`) with both builds.
- Report three numbers per build and thread count: whole-fit wall time; time per fit evaluation
  = wall time / number of evaluations (each is value + 38-parameter gradient); value-only time =
  mean of ten `EvalLLH` calls (`bench_value_only.slurm`). Per evaluation is the fair comparison;
  the whole-fit ratio also contains the different number of evaluations (pitfall 13).
- Figures: `bestfit_expectations.py` once per build (each in its own environment), then
  `plot_accuracy.py` (`PYTHONNOUSERSITE=1`), which writes `timing`, `accuracy_llh`,
  `accuracy_bestfit`, `accuracy_expectations`, `accuracy_bin_differences` and
  `accuracy_same_parameters` to `plots/`.

Reference results (Sapphire Rapids, example fit, tolerances 1e-20):

| | original 1 thread | original 8 threads | fast 1 thread | fast 8 threads |
|---|---|---|---|---|
| fit wall time | 57 min 41 s | 7 min 40 s | 1 min 51 s | 15.5 s |
| evaluations | 1071 | 1071 | 902 | 859 |
| per fit evaluation | 3.23 s | 430 ms | 124 ms | 18 ms |
| per value-only evaluation | 1.06 s | 131 ms | 45 ms | 7.2 ms |
| -ln L at best fit | 20145.560918 | 20145.560918 | 20145.561941 | 20145.561188 |

## 7. Pitfalls met in the reference port

1. **Stale object files.** The original Makefile did not list PhysTools headers as
   dependencies of `GollumFit.o`, so header edits did not recompile it. Two versions of the same
   templates then got linked into one library, which crashed or silently ran old code. Fixed in
   `configure`; still run `make clean` after header changes, and grep the build log for
   `src/GollumFit.cpp` to confirm it recompiled.
2. **Dangling captures in pool tasks.** With a persistent thread pool, a task lambda that
   captures a local lambda by reference outlives it and crashes, but only with more than one
   thread. Capture helpers by value.
3. **Single-threaded runs hide threading bugs**, because the fast engine evaluates inline with
   one thread. Always test two or more threads.
4. **`configure` overwrites the tracked `Makefile`** with machine-specific paths. Never commit the
   generated `Makefile`; change `configure` (and the tracked template, if upstream wants it).
5. **A `sed` that rewrites the prefix twice** produced `local_head_head`, so the build found no
   dependencies. Set the prefix line explicitly and grep it.
6. **Missing Diver library** breaks the Python build at HEAD; patch 02 makes it optional
   (`DIVER_DIR`).
7. **Stale eggs** in a copied prefix's `site-packages` (and `easy-install.pth`) import an old
   module. Delete them and print `GollumFitPy.__file__` in every run.
8. **`~/.local` packages shadow the cvmfs stack** (a broken PIL broke matplotlib). Run Python
   with `PYTHONNOUSERSITE=1`.
9. **gdb fails** with the cvmfs `LD_LIBRARY_PATH`: run
   `env -u LD_LIBRARY_PATH gdb -ex "set environment LD_LIBRARY_PATH $LDP" --args python script.py`.
10. **Baseline provenance.** The user's main checkout had an uncommitted change to the weighters,
    so its installed build was not committed HEAD (its fit gave -ln L 20146.04 instead of
    20145.56). Build the baseline from a clean commit in its own worktree and check it against a
    known result.
11. **Data provenance.** The committed `nullexpectation.npz` (814,847 rows) belongs to an older
    FastMC file than the local `compact.fastmc` (334,683 events). Pair them by row count.
12. **Historical timings in the repository are not baselines.** The committed
    `generate_likelihood_example.err` (2 h 50 min) used the older, larger data and older code.
13. **Evaluation counts are path-dependent at tight tolerances.** Rounding-level differences
    change line-search decisions late in the fit; compare time per evaluation, and use "same
    build and thread count gives identical results" as the reproducibility check.
14. **The login node is a shared single core**; time on compute nodes only.
15. **cvmfs can be unreachable** from the login node; check
    `ls /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_8_x86_64/bin/python` first.
16. **Parameters pass through single-precision `FitParameters`** from Python; steps below about
    1e-6 relative are rounded away.
17. **Do not call `EvalLLH` concurrently on one `GollumFit` object** from several Python threads
    (the GIL is released; the pool and caches are shared). Parallelise MCMC over processes.
18. **Tolerance studies must run on the final build**: an earlier grid ran on an intermediate
    build (pitfall 1) and its evaluation counts were not representative.
19. **Points outside a spline's parameter range** give weight zero to the affected events, in
    both builds; bound proposals or fits to the table extents.
20. **Do not append to registered event containers**; re-register after any reallocation.
21. **After `SetData`, call `ConstructLikelihoodProblem` again** (pseudo-experiment loops).

## 8. Report to the user

State the target commit, route (A or B), every deviation from the reference patches, the gate
results with numbers, the benchmark table of section 6 for the target, where the figures and
logs are, and what was not done. Leave all changes uncommitted unless asked otherwise.
