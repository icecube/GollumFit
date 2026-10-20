#!/bin/bash
#SBATCH --job-name=gollum_injgen
#SBATCH --output=injgen_%j.out
#SBATCH --error=injgen_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=00:30:00
#SBATCH --partition=test
# Generates fake data with the original (HEAD) build. Usage: sbatch --export=ALL,SEED=2027,OUT=injection_corr_data.npz[,EXTRA=--independent] bench_injection_gen.slurm
GOLLUMSPACE=/n/home01/awen/meows_home/GollumFit_test
cd $GOLLUMSPACE/GollumFit-perf/examples/fitting_to_null
source $GOLLUMSPACE/setup_head.sh
export PYTHONNOUSERSITE=1
echo "host $(hostname)"
python generate_injection.py --seed ${SEED:-2027} ${EXTRA} --out ${OUT:-injection_data.npz} 2>&1 >/dev/null | grep "^\[inject\]"
test -s ${OUT:-injection_data.npz}
