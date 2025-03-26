from setuptools import setup, Extension
import sys, os, os.path
import numpy as np
import glob

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

    include_dirs = [gollum_build_path+'/include',
                    COOLCVMFSROOT+'/include',
                    prefix+'/include',
                    np.get_include(),
                    '../include/',
                    '../src/']

    libraries = ['python3', 'boost_python{}{}'.format(sys.version_info[0], sys.version_info[1]),
                 'boost_filesystem','boost_iostreams','boost_system','boost_regex',
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
              extra_compile_args=['-O3','-fPIC','-std=c++11', '-fpermissive'],
              depends=[]),
          ]
      )