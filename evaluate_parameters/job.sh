#!/bin/bash
#SBATCH --time=3000
#SBATCH --mem 48000M
#SBARCH --mem-per-cpu 48000M
#SBATCH --nodes=1
#SBATCH --job-name="TEST48gb"
#SBATCH --cpus-per-task=8
#SBATCH --mail-user emunte@idibell.cat
#SBATCH --mail-type BEGIN        # send email when job begins
#SBATCH --mail-type END          # send email when job ends
#SBATCH --mail-type FAIL         # send email if job fails



echo "Run at:"
# Param
script=$1
params=$2
dataset=$3
include_dir=$4
evaluate_par=$5
log=$6

#echo Rscript $script $params $dataset $include_dir > $log 2>&1
srun --mem 48000M -c 8 Rscript $script $params $dataset $include_dir $evaluate_par > $log 2>&1



echo "Finished" 

