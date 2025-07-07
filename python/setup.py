from setuptools import setup, Extension
import sys, os, os.path
import numpy as np
import glob
import pybind11
import subprocess

if sys.platform == 'win32' or sys.platform == 'win64':
    print('Windows is not a supported platform.')
    quit()

else:
    try:
        gollum_build_path = os.environ['GOLLUMBUILDPATH']
    except KeyError:
        gollum_build_path = '/usr/local/'
    try:
        COOLCVMFSROOT=os.environ['SROOT']
    except KeyError:
        COOLCVMFSROOT= '/usr/local/'
    try:
        prefix=os.environ['PREFIX']
    except:
        prefix = '/usr/local/'

    # add some paths for common locations of libraries and includes
    include_dirs = [gollum_build_path+'/include',
                    COOLCVMFSROOT+'/include',
                    prefix+'/include',
                    np.get_include(),
                    '../include/',
                    '../src/',
                    pybind11.get_include()]

    libraries = ['boost_filesystem','boost_iostreams','boost_system','boost_regex',
                 'LeptonWeighter', 'photospline',
                 'SQuIDS','nuSQuIDS',
                 'gsl','gslcblas','m','z',
                 'hdf5','hdf5_hl','PhysTools','cfitsio']

    architecture = os.uname().machine
    library_dirs = [gollum_build_path+'/lib/python{}.{}/site-packages'.format(sys.version_info[0], sys.version_info[1]),
                    gollum_build_path+'/lib',
                    gollum_build_path+'/lib64',
                    COOLCVMFSROOT+'/lib',
                    COOLCVMFSROOT+'/lib64',
                    prefix+'/lib',
                    prefix+'/lib64',
                    f"/usr/lib/{architecture}-linux-gnu/hdf5/serial"]
    
    def pkgconfig_flags(package, flag):
        try:
            output = subprocess.check_output(['pkg-config', flag, package], encoding='utf-8')
            return output.strip().split()
        except Exception:
            return []
    
    # Add GSL include and library directories if available
    gsl_include_dirs = [flag[2:] for flag in pkgconfig_flags('gsl', '--cflags') if flag.startswith('-I')]
    gsl_library_dirs = [flag[2:] for flag in pkgconfig_flags('gsl', '--libs-only-L') if flag.startswith('-L')]
    gsl_libraries = [flag[2:] for flag in pkgconfig_flags('gsl', '--libs-only-l') if flag.startswith('-l')]
    include_dirs += gsl_include_dirs
    library_dirs += gsl_library_dirs
    libraries += gsl_libraries

    # Add HDF5 include and library directories if available
    hdf5_include_dirs = [flag[2:] for flag in pkgconfig_flags('hdf5', '--cflags') if flag.startswith('-I')]
    hdf5_library_dirs = [flag[2:] for flag in pkgconfig_flags('hdf5', '--libs-only-L') if flag.startswith('-L')]
    hdf5_libraries = [flag[2:] for flag in pkgconfig_flags('hdf5', '--libs-only-l') if flag.startswith('-l')]
    include_dirs += hdf5_include_dirs
    library_dirs += hdf5_library_dirs
    libraries += hdf5_libraries

    # Add cfitsio include and library directories if available
    cfitsio_include_dirs = [flag[2:] for flag in pkgconfig_flags('cfitsio', '--cflags') if flag.startswith('-I')]
    cfitsio_library_dirs = [flag[2:] for flag in pkgconfig_flags('cfitsio', '--libs-only-L') if flag.startswith('-L')]
    cfitsio_libraries = [flag[2:] for flag in pkgconfig_flags('cfitsio', '--libs-only-l') if flag.startswith('-l')]
    include_dirs += cfitsio_include_dirs
    library_dirs += cfitsio_library_dirs
    libraries += cfitsio_libraries

files = ['GollumFitPy.cpp']
gollum_space_path = os.environ['GOLLUMSPACE']
extra_objs = glob.glob(gollum_space_path+'/lib/*.o')

setup(name='GollumFitPy',
      ext_modules = [
          Extension('GollumFitPy', files,
              library_dirs=library_dirs,
              libraries=libraries,
              include_dirs=include_dirs,
              extra_objects=extra_objs,
              extra_compile_args=['-v', '-O3','-fPIC','-std=c++14', '-fpermissive'],
              language='c++',
              depends=[]),
          ]
      )