#!/bin/bash
#SBATCH --job-name=gollum_threads
#SBATCH --output=threads_%j.out
#SBATCH --error=threads_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12000
#SBATCH --time=02:00:00
#SBATCH --partition=arguelles_delgado
# Thread-scaling measurement of the likelihood evaluation on one compute node.
source ../setup.sh
echo "host $(hostname) module $(python -c 'import GollumFitPy; print(GollumFitPy.__file__)')"
for t in 1 2 4 8 16; do
  python bench_llh.py --threads $t --out bench_threads_${t}.npz 2>&1 | grep "^\[bench\]"
done
python bench_compare.py bench_reference.npz bench_threads_1.npz
