#!/bin/bash
#SBATCH --job-name=gollum_valueonly
#SBATCH --output=valueonly_%j.out
#SBATCH --error=valueonly_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --time=01:00:00
#SBATCH --partition=test
# Value-only likelihood evaluation time (EvalLLH, no gradient) for the original (HEAD) and fast builds,
# 1 and 8 threads, 10 repetitions each at the nominal point. Writes timing_value_<build>_t<threads>.npz.
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
export PYTHONNOUSERSITE=1
echo "host $(hostname) $(lscpu | grep "Model name" | sed "s/ \+/ /g")"
for build in head perf; do
  ( source $GOLLUMSPACE/setup_${build}.sh >/dev/null 2>&1
    tag=$([ $build = perf ] && echo fast || echo head)
    for t in 1 8; do
      python bench_llh.py --no-fit --threads $t --repeat 10 --out timing_value_${tag}_t${t}.npz 2>&1 >/dev/null | grep -E "EvalLLH|import"
    done )
done
