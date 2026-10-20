#!/bin/bash
#SBATCH --job-name=gollum_inject
#SBATCH --output=inject_%j.out
#SBATCH --error=inject_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=06:00:00
#SBATCH --partition=test
# Injection-recovery fit, single thread. Usage:
#   sbatch --export=ALL,BUILD=head,PRIORS=truth,DATA=injection_data.npz bench_injection.slurm   (BUILD: head | new; PRIORS: truth | nominal)
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
if [ "${BUILD:-new}" = "head" ]; then source $GOLLUMSPACE/setup_head.sh; else source $GOLLUMSPACE/setup_perf.sh; fi
export PYTHONNOUSERSITE=1
echo "build ${BUILD:-new} priors ${PRIORS} data ${DATA:-injection_data.npz} host $(hostname) $(lscpu | grep "Model name" | sed "s/ \+/ /g")"
/usr/bin/time -f "Elapsed %E, maxRSS %M kB" python fit_injection.py --data ${DATA:-injection_data.npz} --priors ${PRIORS} --threads 1 2>&1 \
  | grep -vE "^Checking file|^Loading|^reset_|^Constructing|^Making|validity|^Remaking|^-?[0-9]|^\[|^Warning"
