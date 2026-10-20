#!/bin/bash
#SBATCH --job-name=gollum_memory
#SBATCH --output=memory_%j.out
#SBATCH --error=memory_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --time=04:00:00
#SBATCH --partition=test
# Memory study: both builds, 1x / 3x / 10x the example's simulated events, 1 and 8 threads.
# Each run is a fresh Python process (mem_probe.py). Results: MEMJSON lines in memory_<job>.out.
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
export PYTHONNOUSERSITE=1
echo "host $(hostname) $(lscpu | grep "Model name" | sed "s/ \+/ /g")"
for build in head perf; do
  ( source $GOLLUMSPACE/setup_${build}.sh >/dev/null 2>&1
    tag=$([ $build = perf ] && echo fast || echo original)
    for k in 1 3 10; do
      f=$([ $k = 1 ] && echo $GOLLUMSPACE/GollumFit-perf/examples/FastMC/compact.fastmc || echo $GOLLUMSPACE/tmp_mem/compact_x$k.fastmc)
      for t in 1 8; do
        python mem_probe.py --compact $f --threads $t --label "$tag x$k t$t" 2>/dev/null | grep "^MEMJSON"
      done
    done )
done
