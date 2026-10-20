#!/bin/bash
#SBATCH --job-name=gollum_fitcmp
#SBATCH --output=fitcmp_%j.out
#SBATCH --error=fitcmp_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --time=08:00:00
#SBATCH --partition=arguelles_delgado
# Runs examples/fitting_to_null/generate_likelihood.py (seed 100, 38 free parameters) with either
# the original build (BUILD=old: GollumFit_test/local, main checkout) or this branch (BUILD=new:
# GollumFit_test/local_perf). Usage:
#   sbatch --export=ALL,BUILD=old,GOLLUM_THREADS=1,GOLLUM_CHANGE_TOL=1e-20,GOLLUM_GRAD_TOL=1e-20 bench_fit_compare.slurm
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
case "${BUILD:-new}" in old) source $GOLLUMSPACE/setup.sh;; head) source $GOLLUMSPACE/setup_head.sh;; *) source $GOLLUMSPACE/setup_perf.sh;; esac
export GOLLUM_THREADS=${GOLLUM_THREADS:-1} GOLLUM_CHANGE_TOL=${GOLLUM_CHANGE_TOL:-1e-20} GOLLUM_GRAD_TOL=${GOLLUM_GRAD_TOL:-1e-20}
echo "build ${BUILD:-new} host $(hostname) threads $GOLLUM_THREADS change_tol $GOLLUM_CHANGE_TOL grad_tol $GOLLUM_GRAD_TOL"
echo "module $(python -c 'import GollumFitPy; print(GollumFitPy.__file__)')"
/usr/bin/time -v python generate_likelihood.py 2> >(grep -E "Elapsed|Maximum resident|User time|signal" >&2) | grep -vE "^Checking file|^Loading|^reset_|^Constructing|^Making|validity|^Remaking|^-[0-9]|^\[" 
