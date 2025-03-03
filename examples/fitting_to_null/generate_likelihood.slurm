#!/bin/bash
#SBATCH --job-name=generate_likelihood_example          
#SBATCH --output=generate_likelihood_example.out       
#SBATCH --error=generate_likelihood_example.err         
#SBATCH --ntasks=1                        
#SBATCH --cpus-per-task=4                 
#SBATCH --mem=10000                       
#SBATCH --time=20:00:00                   
#SBATCH --partition=arguelles_delgado     

source ../../../../setup.sh 

# Run your Python script
srun time python ./generate_likelihood.py