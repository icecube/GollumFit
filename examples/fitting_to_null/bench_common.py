"""Shared configuration for the GollumFit benchmark/regression scripts.

Mirrors generate_likelihood.py (binning, priors, bounds, spline paths) so that
the benchmark exercises exactly the production fit configuration.
"""
import os
import numpy as np

# name: (prior_type, center, width, lower, upper)   -- same table as generate_likelihood.py
SYST = {
    "convNorm":                  ("Gaussian",   1.,   0.2,   0.1,    3.),
    "zenithCorrection":          ("Gaussian",   0.,   1.,   -3.,     3.),
    "kaonLosses":                ("Gaussian",   0.,   1.,   -3.,     3.),
    "hadronicHEkp":              ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicHEkm":              ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE1pip":           ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE1pim":           ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE3kp":            ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE3km":            ("Gaussian",   0.,   1.,   -1.5,    2.),
    "hadronicVHE3pip":           ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE3pim":           ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE3p":             ("Gaussian",   0.,   1.,   -2.,     2.),
    "hadronicVHE3n":             ("Gaussian",   0.,   1.,   -2.,     2.),
    "cosmicRay1":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "cosmicRay2":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "cosmicRay3":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "cosmicRay4":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "cosmicRay5":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "cosmicRay6":                ("Gaussian",   0.,   1.,   -4.,     4.),
    "icegrad0":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad1":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad2":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad3":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad4":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad5":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad6":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad7":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "icegrad8":                  ("Gaussian",   0.,   1.,   -3.,     3.),
    "domEfficiency":             ("Gaussian",   1.27, 0.123, 1.234,  1.346),
    "holeiceForward":            ("Gaussian",  -1.,  10.,   -5.35,   1.85),
    "astroNorm":                 ("Gaussian",   4.72/6., 0.36, 0.,   3.),
    "astroDeltaGamma":           ("Gaussian",   0.,   0.36, -2.,     2.),
    "astroDeltaGammaSec":        ("Gaussian",   0.,   0.36, -2.,     2.),
    "nuxs":                      ("Gaussian",   1.,   0.1,   0.824,  1.176),
    "nubarxs":                   ("Gaussian",   1.,   0.1,   0.824,  1.176),
    "astroPivot":                ("Uniform",    5.,   1.,    4.,     6.),
    "promptNorm":                ("Gaussian",   1.,   1.,    0.,     3.),
    "NeutrinoAntineutrinoRatio": ("Gaussian",   1.,   1.,    0.,     2.),
}
NAMES = list(SYST.keys())

# Parameter points used for the numerical regression check. Each dict is an
# override on top of the nominal (prior-center) point; together they exercise
# every weighter factor (normalisations, gradients, splines, tilt, nu/nubar).
TEST_POINTS = [
    {},
    {"domEfficiency": 1.30},
    {"holeiceForward": 0.5},
    {"nuxs": 1.1, "nubarxs": 0.9},
    {"astroDeltaGamma": 0.3, "astroDeltaGammaSec": -0.2, "astroPivot": 5.3,
     "convNorm": 1.2, "promptNorm": 1.5, "astroNorm": 1.1,
     "icegrad0": 0.5, "icegrad5": -0.7, "hadronicVHE3kp": 0.5, "cosmicRay2": -0.4,
     "kaonLosses": 0.3, "zenithCorrection": 0.2, "NeutrinoAntineutrinoRatio": 1.2,
     "domEfficiency": 1.25, "holeiceForward": -2.0, "nuxs": 0.95, "nubarxs": 1.05},
]


def repo_root():
    """GollumFit root: $GOLLUMFIT_ROOT or two levels above this file."""
    return os.environ.get("GOLLUMFIT_ROOT",
                          os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))


def nominal_params(gf, overrides=None):
    fp = gf.FitParameters()
    for name, (_, center, *_rest) in SYST.items():
        setattr(fp, name, center)
    for k, v in (overrides or {}).items():
        setattr(fp, k, v)
    return fp


def params_to_array(fp):
    return np.array([getattr(fp, n) for n in NAMES], dtype=float)


def make_fit_config(gf, free=()):
    """Flags (True = fixed), bounds and priors as in generate_likelihood.py."""
    flag = gf.FitParametersFlag(True)
    bound = gf.FitParametersBound()
    priors = gf.Priors()
    for name, (ptype, center, width, lo, hi) in SYST.items():
        setattr(bound, name + "Min", lo)
        setattr(bound, name + "Max", hi)
        if ptype == "Gaussian":
            setattr(priors, name + "Center", center)
            setattr(priors, name + "Width", width)
        else:
            setattr(priors, name + "Min", lo)
            setattr(priors, name + "Max", hi)
    for name in free:
        setattr(flag, name, False)
    root = repo_root()
    ice = np.load(os.path.join(root, "resources/correlation_matrices/icegrad_correlations.npy"))
    flux = np.load(os.path.join(root, "resources/correlation_matrices/flux_correlations_new_ddmnodeis.npy"))
    for (i, j), v in np.ndenumerate(ice):
        priors.SetIceGradientsCorr(i, j, v)
    for (i, j), v in np.ndenumerate(flux):
        priors.SetFluxCorr(i, j, v)
    return flag, bound, priors


def make_gollumfit(gf, threads=1, change_tol=1e-20, grad_tol=1e-20, compact=None):
    root = repo_root()
    dp = gf.DataPaths()
    dp.domeff_spline_path = os.path.join(root, "resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99")
    dp.holeice_spline_path = os.path.join(root, "resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99")
    dp.attenuation_spline_path = os.path.join(root, "resources/Splines/AttenuationSplines/new_ddmnodeis")
    dp.compact_file_path = compact or os.path.join(root, "examples/FastMC/compact.fastmc")
    edges = np.logspace(np.log10(300), np.log10(1e5), 25)
    sp = gf.SteeringParams()
    sp.minFitEnergy = edges[0]
    sp.maxFitEnergy = edges[-1]
    sp.logEbinEdge = np.log10(edges[0])
    sp.logEbinWidth = np.log10(edges[1]) - np.log10(edges[0])
    sp.minCosth = -1.0
    sp.maxCosth = 0.0
    sp.cosThbinEdge = 0.0
    sp.cosThbinWidth = 0.05
    sp.selectionStart = 0.99
    sp.evalThreads = threads
    sp.change_tol = change_tol
    sp.grad_tol = grad_tol
    sp.uncertaintyModSigmaOverMu = 0.0
    return gf.GollumFit(dp, sp)


def load_null_data(root=None):
    root = root or repo_root()
    return np.load(os.path.join(root, "examples/fitting_to_null/nullexpectation.npz"))["realization"]


def parse_fits(path):
    """Fits printed by generate_likelihood.py in a log file, in order of appearance.

    Each entry is a dict with keys params (name -> value), llh, neval and wall (MinLLH wall time [s]).
    """
    import re
    fits, cur = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith("MinLLH wall time:"):
            cur = {"params": {}, "wall": float(line.split()[3])}
            fits.append(cur)
            continue
        if cur is None:
            continue
        m = re.match(r"^([A-Za-z0-9]+) (-?[0-9.eE+-]+)$", line)
        if m and m.group(1) in SYST:
            cur["params"][m.group(1)] = float(m.group(2))
        elif line.startswith("llh:"):
            cur["llh"] = float(line.split()[1])
        elif line.startswith("nEval:"):
            cur["neval"] = int(line.split()[1])
    return [f for f in fits if len(f["params"]) == len(SYST) and "llh" in f and "neval" in f]
