#!/bin/bash
#SBATCH --job-name=gollum_head
#SBATCH --output=headfit_%j.out
#SBATCH --error=headfit_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --time=06:00:00
#SBATCH --partition=arguelles_delgado
# Builds committed HEAD (620eec7) in GollumFit_test/GollumFit-head into GollumFit_test/local_head with the
# original flags (configure -> make -> make install -> python setup.py install), then runs the example fit
# (perf worktree's generate_likelihood.py, seed 100, tolerances 1e-20) with 1 and with 8 threads.
set -e
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
source $GOLLUMSPACE/setup_head.sh
cd $GOLLUMSPACE/GollumFit-head
echo "host $(hostname); building HEAD $(git rev-parse --short HEAD) into $PREFIX"
./configure --prefix=$PREFIX --with-boost=$BOOST_DIR --with-nusquids=$PREFIX --with-leptonweighter=$PREFIX --with-squids=$PREFIX --with-phystools=$PREFIX > build_head.log 2>&1
grep -n "^CXXFLAGS +=" Makefile
make -j4 >> build_head.log 2>&1 && make install >> build_head.log 2>&1
cd python && rm -rf build dist *.egg-info && python setup.py install --prefix=$PREFIX > ../build_py_head.log 2>&1
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
echo "module $(python -c 'import GollumFitPy; print(GollumFitPy.__file__)')"
export GOLLUM_CHANGE_TOL=1e-20 GOLLUM_GRAD_TOL=1e-20
for t in 8 1; do
  echo "=== HEAD build, threads $t"
  GOLLUM_THREADS=$t /usr/bin/time -f "Elapsed %E, maxRSS %M kB" python generate_likelihood.py 2>&1 | grep -vE "^Checking file|^Loading|^reset_|^Constructing|^Making|validity|^Remaking|^-?[0-9]|^\["
done
