from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import sys, os, os.path
import numpy as np
import glob
import subprocess
import pybind11

if sys.platform in ("win32", "win64"):
    print("Windows is not a supported platform.")
    sys.exit(1)

gollum_build_path = os.environ.get("GOLLUMBUILDPATH", "/usr/local/")
cvmfs_root = os.environ.get("SROOT", "/usr/local/")
prefix = os.environ.get("PREFIX", "/usr/local/")

include_dirs = [
    gollum_build_path + "/include",
    cvmfs_root + "/include",
    "/usr/local/include",
    prefix + "/include",
    np.get_include(),
    "../include/",
    "../src/",
    pybind11.get_include(),
]

libraries = [
    "boost_filesystem", "boost_iostreams", "boost_system", "boost_regex",
    "LeptonWeighter", "photospline",
    "SQuIDS", "nuSQuIDS",
    "gsl", "gslcblas", "m", "z",
    "hdf5", "hdf5_hl", "PhysTools", "cfitsio", "GollumFit",
    "diver",  # Diver DE optimizer
]

architecture = os.uname().machine
library_dirs = [
    gollum_build_path + f"/lib/python{sys.version_info[0]}.{sys.version_info[1]}/site-packages",
    gollum_build_path + "/lib",
    gollum_build_path + "/lib64",
    cvmfs_root + "/lib",
    cvmfs_root + "/lib64",
    "/usr/local/lib",
    "/usr/local/lib64",
    prefix + "/lib",
    prefix + "/lib64",
    f"/usr/lib/{architecture}-linux-gnu/hdf5/serial",
]

def pkgconfig_flags(package, flag):
    try:
        out = subprocess.check_output(["pkg-config", flag, package], encoding="utf-8")
        return out.strip().split()
    except Exception:
        return []

for pkg in ("gsl", "hdf5", "cfitsio"):
    include_dirs += [f[2:] for f in pkgconfig_flags(pkg, "--cflags") if f.startswith("-I")]
    library_dirs += [f[2:] for f in pkgconfig_flags(pkg, "--libs-only-L") if f.startswith("-L")]
    libraries += [f[2:] for f in pkgconfig_flags(pkg, "--libs-only-l") if f.startswith("-l")]

gollum_space_path = os.environ.get("GOLLUMSPACE", "..")
extra_objs = glob.glob(gollum_space_path + "/lib/*.o")

extra_link_args = []
if sys.platform == "darwin":
    extra_link_args += ["-Wl,-rpath,@loader_path"]
    for d in library_dirs:
        extra_link_args += [f"-Wl,-rpath,{d}"]

ext = Pybind11Extension(
    "GollumFitPy",
    ["GollumFitPy.cpp"],
    library_dirs=library_dirs,
    libraries=libraries,
    include_dirs=include_dirs,
    extra_objects=extra_objs,
    extra_compile_args=["-v", "-O3", "-fPIC", "-std=c++17", "-fpermissive"],
    extra_link_args=extra_link_args,
    language="c++",
)

setup(ext_modules=[ext])