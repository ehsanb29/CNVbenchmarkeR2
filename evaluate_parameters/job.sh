#!/bin/bash
#SBATCH --time=2720
#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --job-name="TEST"
#SBATCH --cpus-per-task=1


echo "Run at:"
# Param
script=$1
params=$2
dataset=$3
include_dir=$4
evaluate_par=$5
log=$6

#echo Rscript $script $params $dataset $include_dir > $log 2>&1
#srun --mem=24G -c 4 
srun --mem=3GB -c 1 Rscript $script $params $dataset $include_dir $evaluate_par > $log 2>&1



echo "Finished" 

