#!/bin/bash
#SBATCH --job-name=gollum_poolmem
#SBATCH --output=poolmem_%j.out
#SBATCH --error=poolmem_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem=40000
#SBATCH --time=01:00:00
#SBATCH --partition=test
# Memory of 16 parallel likelihood workers (fast build, example MC): own object per worker vs forked shared object.
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
source $GOLLUMSPACE/setup_perf.sh
export PYTHONNOUSERSITE=1
echo "host $(hostname) $(lscpu | grep "Model name" | sed "s/ \+/ /g")"
for mode in fork own; do
  python mcmc_fork_memory.py --mode $mode --workers 16 2>/dev/null | grep "^\[pool\]"
done
