# GollumFit likelihood performance: what was changed and why

Branch `perf/fast-likelihood`. This document explains, in order of impact, where the time of a
likelihood fit went, what each change does, how equivalence with the previous code was verified,
and how to benchmark. Numbers refer to the example fit in `examples/fitting_to_null`
(334,683 FastMC meta-events, 24 energy x 20 zenith x 2 topology bins, 38 free nuisance
parameters, L-BFGS-B with autodiff gradients).

## Requirements: matching PhysTools

This branch needs the matching change to PhysTools (the likelihood engine, section 4):

| | branch | commit | based on |
|---|---|---|---|
| GollumFit | `perf/fast-likelihood` | this branch | main 620eec7 |
| PhysTools | `perf/likelihood-engine` | `ced3adf4586900ab2229755c54909fd9aebcedc8` | 8ce44a7 (the commit main uses) |

The `vendor/PhysTools` submodule pointer on this branch is set to `ced3adf`. Until that commit is
pushed to the submodule's remote (`https://github.com/icecube/PhysTools`), `git submodule update`
cannot fetch it; in the meantime check out PhysTools 8ce44a7 and apply
`cookbook/patches/01-phystools-likelihood-engine.patch`, which is the same change. The PhysTools
headers are compiled into GollumFit (PhysTools is header-only for this code), so after switching
PhysTools versions run `make clean` before rebuilding GollumFit.

## 1. Where the time went

A likelihood evaluation loops over every simulated event, computes the event's weight as a
function of the 38 nuisance parameters, sums the weights (and squared weights) per bin, and
feeds the per-bin sums into the SAY likelihood. L-BFGS-B needs the gradient as well, so the
weight was computed with 38-component dual numbers (PhysTools `autodiff::FD<38>`: a value plus
38 derivative slots, 312 bytes). The CPU profile of a full fit (gperftools, symbolised with
pprof) split as:

| Where | Share |
|---|---|
| dual-number arithmetic in the composed weighter chain (`combinedWeighterType::evalImpl`, `FD::operator*`, `operator+`, `autodiff::pow`, memmove of FD temporaries) | ~42 % |
| photospline (`ndsplineeval_core`, `searchcenters`, `bsplvb_simple`, `bspline_deriv_nonzero`) | ~36 % |
| libm (`log10`, `cos`, `pow`, `exp`, `lgamma`) | ~9 % |
| PhysTools per-bin bookkeeping | ~3 % |

The weight of one event is a product of factors:

```
w = convNorm * (conv + prompt + astro) * icegrad                     [enableTotalNorm]
conv    = HoleIce_c * DOMEff_c * Att_c * (1+adu*ad)(1+klu*kl) * F
          F = convFlux + sum_i h_i*had_i + sum_j cr_j*cr_j             (FLT_MAX sentinel if F < 0)
prompt  = promptNorm * HoleIce_p * DOMEff_p * Att_p * promptFlux
astro   = astroNorm  * HoleIce_a * DOMEff_a * Att_a * astroFlux * nn * tilt
          nn   = R (antiparticle) or 2-R (particle)
          tilt = (E/10^pivot)^(-dGamma1 or -dGamma2 depending on E vs 10^pivot)
icegrad = prod_k (1 + g_k * icegrad_k)
```

The expensive factors are the nine spline-based ones: DOM efficiency, hole ice and attenuation
for each of the three flux components. Each is a 3-D photospline table evaluated at
(log10 E, cos zenith, parameter); DOM efficiency and hole ice are applied as
`10^(rate(p) - rate(nominal))`. In the gradient path each factor cost a value and a derivative
evaluation, and the DOM-efficiency and hole-ice factors a third evaluation at the nominal
parameter: 24 full 3-D spline evaluations and 24 `searchcenters` binary searches per event per
likelihood evaluation, plus a `std::map::find` per factor and `log10`/`cos` of event
quantities that never change. On top of that every
multiplication in the chain was a 39-wide dual-number product, with operands passed by value.

Measured before the changes (one thread, login node): a value-only `EvalLLH` took 0.93 s and a
value-plus-gradient evaluation 2.86 s; the committed code needs 1071 evaluations for the
example fit.

## 2. Per-event reduction of the 3-D splines (the big step)

### 2.1 The idea

The part of a spline evaluation that depends on the nuisance parameter is the same for every
event, and the part that depends on the event is the same for every likelihood evaluation. The
original code redid both parts for every event in every evaluation. The new code computes each
part once, where it actually changes, and combines them with three multiply-adds per event.

### 2.2 What a B-spline is

A B-spline is a sum of small, smooth bumps placed along an axis. Each bump `B_k(x)` is a fixed
shape determined only by the knot positions `t`; each has a height `c_k`, the coefficient fitted
when the table was made:

```
f(x) = sum_k c_k * B_k(x)
```

The bumps are piecewise polynomials. Photospline's `order` is the polynomial degree, and the
detector tables here are quadratic (order 2): each bump is made of three quadratic pieces over
three consecutive knot intervals, is zero everywhere else, and the pieces join so that the value
and the slope are continuous. The bumps are defined by the Cox-de Boor recursion, starting from
boxes of degree 0:

```
B[k,0](x) = 1 if t[k] <= x < t[k+1], otherwise 0
B[k,n](x) = (x - t[k]) / (t[k+n] - t[k]) * B[k,n-1](x)
          + (t[k+n+1] - x) / (t[k+n+1] - t[k+1]) * B[k+1,n-1](x)
```

For equally spaced knots one unit apart, the quadratic bump is the familiar bell

```
B(x) = x^2 / 2                on [0, 1)
     = (-2 x^2 + 6 x - 3) / 2  on [1, 2)
     = (3 - x)^2 / 2           on [2, 3),    zero elsewhere
```

Three properties matter here:

- Locality: at any point only three bumps are non-zero (degree + 1), so evaluating `f` needs
  three terms, however many coefficients the table has.
- The bumps are non-negative and add up to one everywhere inside the table range, so `f(x)` is a
  weighted average of the three nearby coefficients, with weights that slide smoothly as `x`
  moves.
- The bump values depend only on `x` and the knots, never on the coefficients.

Evaluating the curve therefore means: find where `x` sits among the knots (a binary search),
compute the three bump values at `x`, and add up three height x bump products.

In three dimensions photospline uses a grid of bumps. Each one is an energy bump times a zenith
bump times a parameter bump, with its own height:

```
f(x0, x1, p) = sum_{i,j,k} c[i,j,k] * B_i(x0) * B_j(x1) * B_k(p)      x0 = log10 E, x1 = cos zenith
```

At any point 3 x 3 x 3 = 27 of these products are non-zero. The coefficient grids are large along
the two event axes and small along the parameter axis:

| table | coefficients (energy x zenith x parameter) | file |
|---|---|---|
| DOM efficiency | 448 x 297 x 8 | 4.3 MB |
| hole ice | 448 x 297 x 7 | 3.7 MB |
| attenuation | 197 x 197 x 27 | 4.2 MB |

### 2.3 Freezing what never changes

In a fit every MC event keeps its energy and zenith; only the nuisance parameter `p` changes from
one evaluation to the next. For a given event `B_i(x0)` and `B_j(x1)` are therefore fixed
numbers, and the sum can be grouped by the parameter bump:

```
f(p) = sum_k [ sum_{i,j} c[i,j,k] * B_i(x0) * B_j(x1) ] * B_k(p)
     = sum_k d[k] * B_k(p)
```

The bracket, `d[k]`, depends only on the event. It is computed once and stored: 8 numbers for
DOM efficiency, 7 for hole ice, 27 for attenuation. What remains is an ordinary one-dimensional
spline in `p` whose heights are the event's own `d[k]`.

A picture: think of the DOM-efficiency table as a stack of 8 maps over the energy-zenith plane,
one per layer along the DOM-efficiency axis. Once per event, read the event's pixel off each map
and write down 8 numbers. For any DOM efficiency after that, blend the 3 neighbouring numbers
with the bump weights. (Strictly, each layer is a map of coefficients rather than of the rate
itself, but the arithmetic is exactly this.)

Nothing is approximated. It is the same sum of products with the brackets placed differently,
the way `a*x + b*x = (a + b)*x`. The only differences are rounding, about one part in 10^6 per
event, which comes from photospline's own single-precision arithmetic.

The parameter bumps `B_k(p)` do not depend on the event, and within one evaluation `p` has a
single value. So the three bump values and their derivatives are computed once per evaluation
and shared by all 334,683 events:

- the event part, `d[k]`, is shared by all evaluations;
- the parameter part, `B_k(p)`, is shared by all events;
- only their product, three multiply-adds (three more for the derivative), is done per event per
  evaluation.

Photospline cannot do this by itself: it is a general evaluator with no memory between calls. It
is handed three coordinates each time and has no way of knowing that two of them are the same as
last time.

### 2.4 What it saves

Per event and per evaluation, the original gradient path needed nine spline factors (DOM
efficiency, hole ice and attenuation for the conventional, prompt and astrophysical fluxes),
each with its value and derivative, and for DOM efficiency and hole ice also the value at the
nominal parameter:

| per event, per evaluation | original | reduced |
|---|---|---|
| full 3-D spline evaluations | 24 | 0 |
| binary searches through the knot grid | 24 | 0 |
| `log10` / `cos` of the event's energy and zenith | about 20 | 0 |
| arithmetic for the spline factors | 24 x 27 products plus the bump values | 9 x 6 multiply-adds |

Memory access is the other half: the original code jumped around 4 MB tables for every event,
which does not fit in the fastest CPU caches; the reduced coefficients are about 500 bytes per
event, stored next to each other. The setup costs about one original evaluation and 170 MB,
paid once per fit against roughly a thousand evaluations.

In the profile of the original fit, spline evaluation was about 36 % of the time and the
38-component dual-number gradient about 42 %. The spline reduction removes the first share and
most of the `log10`/`cos` cost, so on its own it would give roughly a 1.5-2x speed-up (an
estimate from the profile shares, not a separate measurement). The large speed-up comes from
pairing it with the closed-form gradient of section 3, and the pairing is natural: the reduced
spline gives the derivative with respect to `p` in the same three multiply-adds, which is
exactly what the closed-form gradient needs.

### 2.5 Implementation

`include/ReducedSplines.h` (`ReducedAxisSpline`) implements this with photospline's public API:
it reproduces `splinetable::searchcenters` for the two fixed axes, calls photospline's own
`bsplvb_simple` for the basis values, and contracts the coefficient table over the two fixed
axes (accumulating in double, storing float). At evaluation time the parameter value is the
same for every event, so the basis functions along p (order+1 = 3 non-zero values and their
derivatives, from `bsplvb_simple` / `bspline_deriv_nonzero`) are computed once per table per
evaluation (`ReducedAxisSpline::basis`). Each event then costs three fused multiply-adds for
the value and three for the derivative per factor (`eval`, `evalDeriv`).

`include/SystematicsCache.h` holds the per-event rows: for each of the three flux components
the reduced DOM-efficiency, hole-ice and attenuation coefficients, the log-rate at the nominal
DOM efficiency (1.27) and hole-ice parameter (-1.0, the values hard-coded in the original
weighters), and per-event flags: neutrino/antineutrino (which attenuation scale applies),
particle/antiparticle (nu/nubar balance), whether attenuation is skipped (cos primary zenith
> 0.1, or no table for this particle type), and ln E for the astrophysical tilt. It costs about
500 bytes per event (170 MB for the example) and is built lazily the first time a weighter is
constructed. Events are located by address, which is why the event containers became
`std::vector<Event>` instead of `std::deque<Event>` (`GollumFit::RegisterSimulationWithWeighter`
tells the weighter maker which containers exist).

The edge semantics of the original weighters are reproduced: an event outside a table's
support or with zero nominal rate gets weight 0; a missing attenuation table or a down-going
event gets factor 1; a negative conventional flux returns the FLT_MAX sentinel with zero
gradient; a parameter outside the table range makes the factor 0. Because the rate and the
nominal rate now come from the same reduced coefficients, the DOM-efficiency and hole-ice
corrections are exactly 1 at the nominal point (the committed behaviour; the uncommitted
`cachedDOMEff*` variant on main gave corrections != 1 at nominal for FastMC meta-events).

## 3. Closed-form gradient instead of dual numbers

Each factor of the weight depends on only one or two parameters, so the gradient is sparse and
has a closed form, `dw/dp_i = w * d ln f_i / dp_i`. `include/FastWeighter.h`
(`FastWeighterCore::evaluate`) computes the weight once and then writes all 38 partial
derivatives directly:

- normalisations: `dw/dconvNorm = (conv+prompt+astro)*icegrad`; promptNorm and astroNorm
  use the component without its normalisation;
- the 16 hadronic / cosmic-ray gradient parameters: the common factor
  `convNorm * HoleIce_c * DOMEff_c * Att_c * B * icegrad` times the event's cached value
  (zero if the FLT_MAX sentinel is active);
- the atmospheric density and kaon-loss parameters: the same factor with the other
  `(1 + x)` term;
- ice gradients: `S * icegrad_k * prod_{j != k} (1 + g_j c_j)`, computed with prefix and
  suffix products so that nothing is divided by a factor that could be zero;
- DOM efficiency and hole ice: the factor is `10^(r(p) - r0)`, so its derivative is
  `ln10 * factor * r'(p)` with `r'` from the reduced spline;
- the astrophysical tilt `exp(-dGamma (ln E - pivot ln10))` gives the dGamma and pivot
  derivatives in closed form, with the branch chosen by whether E is above the pivot;
- nu/nubar balance: +-1 times the astro component without the balance factor;
- attenuation: the derivative goes into the nuxs slot for neutrinos and the nubarxs slot
  for antineutrinos.

That is about 150 flops and six `exp` calls per event instead of several thousand flops and
seven dual-number `pow` calls. The result is still packed into one `FD<38>` per event
(`FastGradWeighter`), so the summation in PhysTools and the SAY likelihood are unchanged.
`sterile::WeighterMaker::operator()` now returns `FastWeighter` (for `std::vector<double>`)
or `FastGradWeighter` (for `std::vector<FD<38>>`); the original composed chain is kept as
`WeighterMaker::LegacyWeighter<T>()` and is the definition of the weight.

## 4. Likelihood engine (vendored PhysTools), build flags, stopping criteria

In `vendor/PhysTools/PhysTools/likelihood/likelihood.h`:

- a `ThreadPool` was created and joined on every evaluation; it is now persistent
  (`ThreadPool::size()` added);
- per-event weights were pushed into two `std::vector<FD<38>>` per bin and only then
  Kahan-summed; they are now accumulated with the same compensated summation in running form;
- every bin's entry vector was deep-copied about eleven times per evaluation by the histogram
  iterators, by-value returns and casts; observation and simulation bins are now paired once
  (`BinPairCache`, using a new `histogram::binRef`) and grouped into load-balanced tasks
  (longest bins first, four tasks per thread); single-threaded evaluation runs inline;
- the dual number was printed to stdout after every evaluation; removed. The BAD BIN event
  dump (which produced tens of MB when a fit strayed) is behind `setVerboseBadBins(true)`;
- `BFGS_Function` copied both histograms per fit; it now holds a reference;
- `gammaPriorPoissonLikelihood` used a heap vector per bin; stack array.

Build flags: `-O3 -DNDEBUG -march=haswell -fno-math-errno` in `configure`, `Makefile`,
`python/setup.py` and PhysTools' `build/config.mk` (GCC 9.2 has no `x86-64-v3`). `-DNDEBUG`
removes three asserts per event; `-march=haswell` enables AVX2/FMA for the 39-wide loops;
`-ffast-math` is deliberately not used (the code relies on NaN checks and sentinels).

Stopping criteria: the driver's `change_tol` is an absolute change of the log-likelihood
between iterations and `grad_tol` the largest scaled projected-gradient component (scale 1e-3
per parameter, entering only this test). Both defaulted to 1e-20, which means "run until
floating-point noise". Section 6.5 shows that the change criterion is the active one, that 1e-6
reproduces the 1e-20 fit (likelihood within 5e-5 on the final build) with about 7 % fewer
evaluations, and that anything looser stops early along the nearly flat holeiceForward
direction. `SteeringParams` now defaults to 1e-6 / 1e-6.

Other fixes made on the way: `SetData` invalidates the likelihood problem (it held references
to the destroyed data events); `FitResult` counters are initialised; the compact file is read
in one block; the Diver binding is optional (`DIVER_DIR`); the GIL is released around the
heavy Python calls.

## 5. Verification

- **Same code paths, same parameters.** `examples/fitting_to_null/bench_llh.py --no-fit` stores the
  bin expectations and likelihoods at five parameter points that together exercise every factor of
  the weight. The snapshot of the original build (`bench_head.npz`, committed HEAD built from
  scratch, see section 6.1) is bit-identical to the earlier baseline snapshot `bench_reference.npz`.
  Against it the fast build (`bench_final.npz`) agrees to 2.5e-6 relative in every bin. The
  likelihoods agree to 1e-6 at the nominal point and to 1.4e-3 at the farthest test point
  (LLH 23738); one sigma is 0.5. Compare with `bench_compare.py bench_head.npz bench_final.npz`.
- **Every event.** `GollumFit::CheckFastWeighter(fp)` (Python `CheckFastWeighter`,
  `check_fast_weighter.py`) evaluates every event with the fast and the legacy weighter, value and
  all 38 derivatives. Differences are of order 1e-5 per event, the single-precision rounding of
  photospline's own evaluation.
- **Whole fits.** The example fit with both builds, section 6.
- **Tests.** `pytest tests/test_gollumfit.py`: 42/42.

## 6. Results: original build vs fast build

### 6.1 Setup

| | original | fast |
|---|---|---|
| code | committed HEAD 620eec7, worktree `GollumFit_test/GollumFit-head` | branch `perf/fast-likelihood`, worktree `GollumFit_test/GollumFit-perf` |
| install prefix / environment | `local_head` / `setup_head.sh` | `local_perf` / `setup_perf.sh` |
| compiler flags | `-std=c++17 -O3` (as configured on main) | `-std=c++17 -O3 -g -DNDEBUG -march=haswell -fno-math-errno` |
| built by | `cookbook/slurm/build_and_fit_head.slurm` (configure, make, install, `setup.py install`) | the same steps |

The only change to HEAD is that the Python build no longer requires the Diver library
(`python/setup.py`, `python/GollumFitPy.cpp`); the likelihood code is untouched.

Both builds ran the same script, `examples/fitting_to_null/generate_likelihood.py`, with seed 100
(the initial values of all 38 parameters are identical in every log), the same FastMC file
`examples/FastMC/compact.fastmc` (334,683 meta-events) and fake data `nullexpectation.npz`. The
tolerances were the committed example's 1e-20, so both minimisers run until they stop making
progress. All four fits ran on the FASRC `test` partition on one Intel Xeon Platinum 8480CL
(Sapphire Rapids) node, `holy8a24601`; the value-only timings ran on `holy8a24402` of the same type.
The fast single-thread fit overlapped with the original build's compilation and 8-thread fit on
that node; its time matches an earlier run on a Cascade Lake node within the expected
hardware difference. Fit wall time is `MinLLH` only; loading the data adds about 2 s.

### 6.2 Timing

Figure: `examples/fitting_to_null/plots/timing.pdf`, three panels: whole fit, one fit
evaluation, one value-only evaluation, with times written as h/min/s or ms.

| build | threads | fit wall time | evaluations | per fit evaluation | per value-only evaluation | LLH |
|---|---|---|---|---|---|---|
| original | 1 | 3461 s (58 min) | 1071 | 3.23 s | 1.06 s | 20145.560918 |
| original | 8 | 460 s | 1071 | 0.430 s | 0.131 s | 20145.560918 |
| fast | 1 | 111 s | 902 | 0.124 s | 0.045 s | 20145.561941 |
| fast | 8 | 15.5 s | 859 | 0.018 s | 0.0072 s | 20145.561188 |

| speed-up, fast over original | 1 thread | 8 threads |
|---|---|---|
| whole 38-parameter fit | 31x | 30x |
| per fit evaluation (value + gradient) | 26x | 24x |
| per value-only evaluation | 23x | 18x |

How the columns are defined:
- **Per fit evaluation** is the fit wall time divided by the number of evaluations. L-BFGS-B only
  ever requests the value together with the 38-parameter gradient, so every evaluation in a fit
  computes both; the minimiser's own work for 38 parameters is negligible.
- **Per value-only evaluation** is the mean of ten `EvalLLH` calls at the nominal point
  (`bench_value_only.slurm`, `timing_value_<build>_t<threads>.npz`). The gradient makes an
  evaluation 3.0x (original) or 2.7x (fast) as expensive as the value alone.
- **Whole fit** also includes the different number of evaluations (section 6.4), so per
  evaluation is the cleaner measure of the code speed.

Going from 1 to 8 threads, the original build speeds up 7.5x (fit evaluation) and 8.1x (value
only), the fast build 6.8x and 6.3x. At 7-18 ms per evaluation the serial parts of an evaluation
(building the weighter and the spline bases, the prior, handing out tasks) start to show.

### 6.3 Accuracy

| fit | LLH | difference to original, 1 thread | largest parameter shift (prior widths) |
|---|---|---|---|
| original, 1 thread | 20145.560918 | reference | |
| original, 8 threads | 20145.560918 | identical in every digit and parameter | 0 |
| fast, 1 thread | 20145.561941 | +1.0e-3 | astroPivot -0.077, then cosmicRay1 +0.021 |
| fast, 8 threads | 20145.561188 | +2.7e-4 | astroPivot -0.078, then hadronicVHE1pim -0.012 |

Figures (`examples/fitting_to_null/plots/`), all comparing against the original build's
single-thread fit:
- **`accuracy_llh.pdf`**: a table of the four best-fit values of -ln L and their difference to
  the original single-thread fit (0, 0, +1.0e-3, +2.7e-4; one sigma for one parameter is 0.5).
- **`accuracy_bestfit.pdf`**: top, the best-fit value of each of the 38 parameters for all four
  fits, as (fit - prior center) / prior width, with the original single-thread fit as the
  baseline; bottom, each fit minus that baseline, in prior widths. Everything is within 0.025
  prior widths except astroPivot at -0.08. astroPivot has a flat uniform prior and the likelihood
  is nearly flat along it, so a fit that stops at slightly different rounding lands at a
  slightly different point.
- **`accuracy_expectations.pdf`**: expected track spectra in energy and cos zenith at each
  fit's best-fit point, over the fake data; lower panels show the ratio to the original
  single-thread fit. The fast fits agree to better than 1e-4 below 10 TeV and to 7e-4 in the
  highest-energy bins, which contain about 10 events.
- **`accuracy_bin_differences.pdf`**: per-bin relative difference of the best-fit expectation,
  fast vs original, left at 1 thread (median 4e-5, max 1.0e-3) and right at 8 threads
  (median 6e-5, max 1.3e-3), over all 960 bins.
- **`accuracy_same_parameters.pdf`**: the code-only comparison: per-bin relative difference at
  identical parameters, all bins at the five test points (max 2.5e-6).

The likelihood differences are 1e-3 or less on a scale where one sigma is 0.5, and all shifts
are along the flattest directions of the fit, so the two builds give the same physics result.

### 6.4 Why the numbers of evaluations differ

The seeds are identical, and the fits are deterministic: the same build at the same thread
count gives the same likelihood to every digit and the same number of evaluations, even on
different hardware (fast build: 902 evaluations at 1 thread and 859 at 8 threads on both a
Cascade Lake and a Sapphire Rapids node). The counts differ because the arithmetic differs at the
rounding level. The fast build sums the same terms in a different order, its spline factors
differ from the original's at about 1e-6, and it groups bins differently for 1 and 8 threads.
With tolerances of 1e-20 the minimiser runs until it can no longer improve; the last few hundred
steps wander along nearly flat directions (astroPivot, holeiceForward), where a rounding-level
difference can change one line-search decision, after which the paths separate and stop at
different points. The original build sums the bins in the same order for any thread count, which
is why its 1- and 8-thread fits are identical.

### 6.5 Stopping tolerances

`SteeringParams` now defaults to `change_tol = grad_tol = 1e-6` (section 4). On the final build
(Cascade Lake nodes, jobs 51157924, 51157926, 51157928 and 51157929):

| threads | 1e-20: evaluations, LLH | 1e-6: evaluations, LLH | difference |
|---|---|---|---|
| 1 | 902, 20145.5619406 | 844, 20145.5619410 | 6 % fewer evaluations, dLLH 5e-7 |
| 8 | 859, 20145.5611879 | 795, 20145.5612374 | 7 % fewer evaluations, dLLH 5e-5 |

A wider grid was run earlier (`tolstudy_*.log`). **It ran on an intermediate build**: the library
object had not been recompiled after the last engine change, because of the missing header
dependency described in section 7. Its evaluation counts are therefore not those of the final
build (1253 evaluations at 1e-20 instead of 902), but its qualitative result stands. The
absolute change of the likelihood is the active criterion, and anything looser than 1e-6 stops
early along the nearly flat holeiceForward direction (prior width 10):

| change_tol / grad_tol | evaluations | dLLH to 1e-20 | largest parameter shift |
|---|---|---|---|
| 1e-6 / 1e-6 | 983 | +1.1e-3 | holeiceForward 0.04 |
| 1e-6 / 1e-4 | 819 | +1.4e-2 | holeiceForward 0.25 |
| 1e-5 / 1e-6 | 456 | +0.13 | holeiceForward 2.5 |
| 1e-4 / 1e-5 | 225 | +0.25 | holeiceForward 2.8 |
| 1e-3 / 1e-4 | 121 | +0.78 | holeiceForward 2.9 |

Fewer evaluations would need better conditioning of the parameter space rather than looser
tolerances: the per-parameter `scale` only enters the stopping test, so rescaling parameters by
their prior widths inside `MinLLH`, or a preconditioned minimiser, is the next lever.

### 6.6 Thread scaling of the fast build

Cascade Lake node, `bench_threads.slurm`, job 51157936 (`threads_51157936.out`):

| threads | value-only evaluation | fit evaluation (value + gradient) | speed-up |
|---|---|---|---|
| 1 | 0.056 s | 0.142 s | 1.0 |
| 2 | 0.028 s | 0.072 s | 2.0 |
| 4 | 0.014 s | 0.036 s | 3.9 |
| 8 | 0.007 s | 0.018 s | 7.9 |
| 16 | 0.005 s | 0.010 s | 14 |

### 6.7 Memory

The fast build trades memory for speed: the per-event table of reduced spline coefficients
(`SystematicsCache`) is kept for every simulated event. Measured with `mem_probe.py`
(`bench_memory.slurm`, job 52993200, Sapphire Rapids): the resident memory of a likelihood
session at 1x, 3x and 10x the example's 334,683 simulated events; the larger sets repeat the
example's events (`make_scaled_compact.py`). Figure: `plots/memory_scaling.pdf`.

| simulated events | original, peak | fast, peak | extra |
|---|---|---|---|
| 334,683 (the example) | 394 MB | 550 MB | +156 MB |
| 1.00 M | 652 MB | 1139 MB | +487 MB |
| 3.35 M | 1557 MB | 3201 MB | +1644 MB |

Single thread; at 8 threads the numbers are the same within 50 MB (the original is slightly
higher there, from its per-bin gradient buffers). Memory is linear in the number of simulated
events, with the same fixed part of about 260 MB in both builds:

| per simulated event | original | fast |
|---|---|---|
| event record (`Event`, 344 B) and container overhead | 360 B | 344 B |
| histogram references, data, likelihood problem | 24 B | 24 B |
| per-event spline coefficients (`SystematicsCache`: 126 floats + 40 B of flags) | 0 | 544 B |
| total, measured slope | 405 B | 923 B |

The extra memory is set by the number of simulated events alone. It does not grow with the
number of data events, bins, threads, fit parameters or evaluations, and nothing accumulates
during a fit. The table is built on the first evaluation (7 s for 3.35 M events), which is where
the jump in the figure happens. Of the 544 B, 324 B are the attenuation coefficients (27 per flux
component), 168 B DOM efficiency and 84 B hole ice; for an exact method this is close to the
minimum (half-precision storage would cost accuracy).

Extrapolated with these slopes:

| use | simulated events | original | fast |
|---|---|---|---|
| example FastMC | 0.33 M | 0.4 GB | 0.55 GB |
| y-tau full-MC fit | 1.03 M | 0.7 GB | 1.2 GB |
| full BDT MC without FastMC | 13.1 M | 5.3 GB | 11.8 GB |

Where it could bite:
- **Very large MC fitted directly** (without FastMC): 13 M events need about 12 GB instead of 5 GB.
- **Loading raw MC from HDF5** (`LoadMC`): the event container is now a contiguous
  `std::vector` filled one event at a time; when it grows it briefly holds the old and the new
  copy, so loading can transiently need up to twice the 344 B per event. Not measured here: the
  compact-file path reads all events in one step and has no such peak.
- **`ConstructFastMode` followed by fitting in the same process**: both the full MC and the
  meta-events are registered with the weighter, so the table is built for the full MC as well,
  although only the meta-events are fitted. The usual workflow (write the compact file, fit it
  in a new job) does not hit this.
- **Several fits in parallel processes** (e.g. MCMC walkers): each process has its own copy,
  e.g. 8 x 1.2 GB for the y-tau MC, unless the processes are forked from one prepared object
  (below).

**Parallel MCMC.** Memory scales with the number of worker processes, not with the number of
walkers: an ensemble sampler keeps each walker as 38 numbers and a pool of P processes evaluates
them in turn. Worker processes forked from one prepared object share its memory, because the
event records and the per-event table are only read during evaluation. Measured with
`mcmc_fork_memory.py` (`bench_mcmc_memory.slurm`, job 53955528), 16 workers on the example MC,
fast build; "physical memory" is the sum of PSS over all processes:

| pool of 16 workers | physical memory | per extra worker |
|---|---|---|
| every worker builds its own object | 8.3 GB | 517 MB |
| workers forked from one prepared object | 0.59 GB | 3 MB |

Every worker reproduced the parent's likelihood exactly. Tools that add up per-process RSS
(`top`, some accounting) report about 9 GB for the forked pool, because they count the shared
pages once per process.

Recipe (emcee or any sampler that takes a `multiprocessing` pool):

```
import multiprocessing as mp
gollumfit = build_gollumfit(evalThreads=1)   # module-level object; one evaluation thread is required
gollumfit.EvalLLH(start_params, True)        # warm-up: builds the per-event table before forking

def log_prob(x):                              # uses the module-level object; it cannot be pickled
    return -gollumfit.EvalLLH(to_fit_parameters(x), True)

with mp.get_context("fork").Pool(n_cores) as pool:
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob, pool=pool)
```

Three conditions:
- **Fork.** Use the `fork` start method explicitly. It is the default on Linux up to Python 3.13;
  `spawn` and `forkserver` would build a new object in every worker.
- **One evaluation thread.** Threads do not survive a fork, so a pool of evaluation threads
  created in the parent would leave the workers waiting for threads that no longer exist.
- **Warm-up before the pool.** Without the warm-up evaluation, every worker would build its own
  per-event table and the sharing would be lost.

Possible mitigations, not implemented: build the table only for the container the simulation
histogram uses (removes the third case); reserve the container before `LoadMC` fills it, or fill
it once at the end (removes the loading peak); a steering switch to use `LegacyWeighter` in
memory-limited jobs.

### 6.8 Injection-recovery test

Fake data (Asimov, made with the original build) at a random point in the 38-parameter space,
fitted by both builds (single thread, tolerances 1e-20, the example's seeded starting point), with
the Gaussian priors centred either on the injected values or at nominal. Two injections:
"correlated", drawn from the priors as the fit uses them, including the flux and ice-gradient
correlations (`generate_injection.py`, seed 2027); "independent", each parameter drawn on its own
(`--independent`, seed 2026). The flux correlation matrix is nearly degenerate (smallest
eigenvalue 7e-4), so the independent point sits about 60 sigma from the prior; with nominal priors
that case is a stress test rather than a realistic one. Scripts `fit_injection.py`,
`bench_injection.slurm`, `bench_injection_gen.slurm`; summary `plot_injection.py`. Figures:
`plots/injection_summary.pdf`, `plots/injection_recovery_injection_corr_data.pdf`,
`plots/injection_recovery_injection_data.pdf`.

| injection, priors at | -ln L at truth | original: best fit, evaluations, time | fast: best fit, evaluations, time |
|---|---|---|---|
| correlated, truth | 20307.9889 | 20307.9616, 1222, 1 h 7 min | 20307.9724, 1119, 2 min 10 s |
| correlated, nominal | 20325.3772 | 20313.1093, 527, 29 min | 20313.0945, 686, 1 min 23 s |
| independent, truth | 20199.8879 | 20199.8237, 1453, 1 h 18 min | 20199.8237, 1472, 3 min 8 s |
| independent, nominal | 22205.6765 | 20203.7508, 1194, 1 h 5 min | 20203.7516, 1292, 2 min 40 s |

- **Agreement between the builds**: the best-fit -ln L differ by at most 0.015 (one sigma is 0.5),
  with neither build systematically lower, and the best-fit parameters by at most 0.08 prior
  widths (cosmicRay1).
- **Recovery**: with priors centred at the truth, both builds find a point slightly better than
  the truth itself and recover every parameter within 0.11 prior widths, except astroPivot in the
  correlated injection (0.66 for both builds: the injected spectral indices below and above the
  break are nearly equal, 0.104 and 0.096, so the break energy is barely constrained) and the
  neutrino/antineutrino ratio in the independent one (0.27 for both).
- **With nominal priors** the fits land between the truth and the prior centres, the same way for
  both builds.
- The `succeeded` flag differs in one case (fast False, original True for the correlated
  injection with nominal priors). With tolerances of 1e-20, L-BFGS-B normally stops when its line
  search can no longer make progress and reports that as a failure; the flag follows the path,
  not the quality of the minimum.

## 7. Files

All paths below are relative to `/n/home01/awen/meows_home/GollumFit_test/`.

| file | what it is |
|---|---|
| `GollumFit-perf/examples/fitting_to_null/plots/timing.pdf` | timing figure, section 6.2 |
| `GollumFit-perf/examples/fitting_to_null/plots/accuracy_llh.pdf` | table of the four best-fit likelihood values, section 6.3 |
| `GollumFit-perf/examples/fitting_to_null/plots/accuracy_bestfit.pdf` | best-fit parameters and differences to the baseline, section 6.3 |
| `GollumFit-perf/examples/fitting_to_null/plots/accuracy_expectations.pdf` | best-fit track spectra and ratios, section 6.3 |
| `GollumFit-perf/examples/fitting_to_null/plots/accuracy_bin_differences.pdf` | per-bin best-fit differences at 1 and 8 threads, section 6.3 |
| `GollumFit-perf/examples/fitting_to_null/plots/accuracy_same_parameters.pdf` | per-bin differences at identical parameters (code only), section 6.3 |
| `GollumFit-perf/examples/fitting_to_null/plots/memory_scaling.pdf` | memory vs number of simulated events, section 6.7 |
| `GollumFit-perf/examples/fitting_to_null/memory_52993200.out` | memory study results (MEMJSON lines, read by `plot_memory.py`) |
| `GollumFit-perf/examples/fitting_to_null/poolmem_53955528.out` | parallel-worker memory, own object vs forked (section 6.7) |
| `GollumFit-perf/examples/fitting_to_null/plots/injection_summary.pdf` | injection-recovery results, section 6.8 |
| `GollumFit-perf/examples/fitting_to_null/inject_*.out`, `injection_data.npz`, `injection_corr_data.npz` | injection fits and their fake data (section 6.8) |
| `GollumFit-head/headfit_51361096.out` (`.err`) | original build: build log line, then the 8-thread fit, then the 1-thread fit |
| `GollumFit-perf/examples/fitting_to_null/fitcmp_51327941.out` (`.err`) | fast build, 1-thread fit |
| `GollumFit-perf/examples/fitting_to_null/fitcmp_51327946.out` (`.err`) | fast build, 8-thread fit |
| `GollumFit-perf/examples/fitting_to_null/valueonly_51362702.out` | value-only timings, both builds |
| `GollumFit-perf/examples/fitting_to_null/timing_value_{head,fast}_t{1,8}.npz` | value-only timings (`t_eval`: ten repetitions) |
| `GollumFit-perf/examples/fitting_to_null/bestfit_{head,fast}.npz` | best-fit parameters, expectations and recomputed LLH of the four fits |
| `GollumFit-perf/examples/fitting_to_null/bench_head.npz`, `bench_final.npz` | fixed-parameter snapshots of the original and fast builds |
| `GollumFit-perf/examples/fitting_to_null/tolstudy_*.log` | earlier tolerance grid (intermediate build, section 6.5) |
| `GollumFit-perf/examples/fitting_to_null/threads_51157936.out` | thread scaling, section 6.6 |
| `GollumFit-head/build_head.log`, `build_py_head.log` | build logs of the original build |

Superseded and kept only for the record: `fitcmp_51157924` to `fitcmp_51157934` on Cascade Lake
nodes, whose "old" build was the main checkout's working tree with the uncommitted
`cachedDOMEff*` change rather than committed HEAD, and `bench_reference.npz`, identical to
`bench_head.npz`. The 2 h 50 min in the committed `generate_likelihood_example.err` is not a
comparable baseline either: it used the older example data (814,847 FastMC events, 2.4 times
the current set) with older code.

A build caveat: before 2026-10-07, `GollumFit.o` did not depend on the PhysTools headers in
the Makefile, so editing those headers did not recompile it. This is fixed in `configure` and
`Makefile`; after pulling header changes made elsewhere, run `make clean` once.

## 8. Reproduction

```
cd /n/home01/awen/meows_home/GollumFit_test
sbatch -p test cookbook/slurm/build_and_fit_head.slurm         # original build + its 8- and 1-thread fits
cd GollumFit-perf/examples/fitting_to_null
sbatch -p test --export=ALL,BUILD=new,GOLLUM_THREADS=1 bench_fit_compare.slurm   # fast build fits
sbatch -p test --export=ALL,BUILD=new,GOLLUM_THREADS=8 bench_fit_compare.slurm   # (BUILD=head also works)
sbatch bench_value_only.slurm                                   # value-only timings, both builds

export PYTHONNOUSERSITE=1                                       # a broken PIL in ~/.local breaks matplotlib
source ../../../setup_head.sh
python bench_llh.py --no-fit --out bench_head.npz               # fixed-parameter snapshot, original build
python bestfit_expectations.py --out bestfit_head.npz head8=../../../GollumFit-head/headfit_<job>.out:0 head1=../../../GollumFit-head/headfit_<job>.out:1
source ../../../setup_perf.sh
python bench_llh.py --no-fit --out bench_final.npz              # fixed-parameter snapshot, fast build
python bestfit_expectations.py --out bestfit_fast.npz fast1=fitcmp_<job>.out fast8=fitcmp_<job>.out
python bench_compare.py bench_head.npz bench_final.npz         # regression gate
python check_fast_weighter.py                                   # fast vs legacy weighter, every event
python plot_accuracy.py                                         # all six figures (.pdf and .png)
```

## 9. Porting to `inelasticity_dev` and other installs

`cookbook/COOKBOOK.md` is the step-by-step procedure, with the changes as three patches in
`cookbook/patches/` (PhysTools engine; GollumFit core; examples, benchmarks and docs). The
patches apply cleanly to main (620eec7) and PhysTools 8ce44a7. On `inelasticity_dev`, PhysTools
is the same commit and `analysisWeighting.h` differs by one include, but `Event.h`,
`GollumParameters.h`, `GollumFit.h`, `src/GollumFit.cpp`, `compactIO`, `configure` and the
Python bindings differ, so the core patch has to be applied by hand there (cookbook section
4.4). The side table costs ~500 B per event (~0.5 GB for a 1 M-event full-MC fit).
