from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import sys, os, os.path
import numpy as np
import glob
import subprocess
import pybind11
from pathlib import Path

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

def env_flag(name):
    value = os.environ.get(name)
    if value is None:
        return None
    value = value.strip().lower()
    if value in ("1", "true", "yes", "on"):
        return True
    if value in ("0", "false", "no", "off"):
        return False
    raise RuntimeError(f"{name} must be one of 1/0, true/false, yes/no, or on/off")

def find_diver(diver_dir=None):
    roots = []
    if diver_dir:
        roots.append(Path(diver_dir).expanduser())
    roots += [Path(prefix), Path(gollum_build_path), Path("/usr/local")]

    seen = set()
    for root in roots:
        root = root.resolve()
        if root in seen:
            continue
        seen.add(root)

        include_dir = root / "include"
        if not (include_dir / "diver.hpp").is_file():
            continue

        for lib_dir in (root / "lib", root / "lib64", root):
            if any((lib_dir / name).is_file() for name in ("libdiver.so", "libdiver.a", "libdiver.dylib")):
                return include_dir, lib_dir

    return None, None

with_diver = env_flag("GOLLUMFIT_WITH_DIVER")
diver_dir = os.environ.get("DIVER_DIR")
define_macros = []
diver_include, diver_lib = find_diver(diver_dir)
diver_requested = with_diver is True or bool(diver_dir)
diver_found = diver_include is not None and diver_lib is not None

if with_diver is False:
    print("Building GollumFitPy without Diver support because GOLLUMFIT_WITH_DIVER=0.")
elif diver_found:
    include_dirs.append(str(diver_include))
    library_dirs.append(str(diver_lib))
    libraries.append("diver")
    define_macros.append(("GOLLUMFIT_HAVE_DIVER", "1"))
    print(f"Building GollumFitPy with Diver support from {diver_include.parent}")
elif diver_requested:
    hint = "Set DIVER_DIR to a Diver checkout containing include/diver.hpp and lib/libdiver.so, or install those files under PREFIX/include and PREFIX/lib."
    raise RuntimeError("Diver support was requested but Diver was not found. " + hint)
else:
    print("Building GollumFitPy without Diver support. Install Diver into PREFIX or set DIVER_DIR=/path/to/Diver to enable it.")

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
    define_macros=define_macros,
    extra_compile_args=["-v", "-O3", "-fPIC", "-std=c++17", "-fpermissive"],
    extra_link_args=extra_link_args,
    language="c++",
)

setup(
    name='GollumFitPy',
    version='1.2.0',
    ext_modules=[ext]
)
