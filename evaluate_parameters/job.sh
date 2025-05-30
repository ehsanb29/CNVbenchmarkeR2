#$ -S /bin/bash
#$ -e /myPath/output.txt
#$ -o /myPath/output.txt

# Determine number of threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Create temp dir for GATK
mkdir /tmp/tmpGATK
tmp_dir=$(mktemp -p /tmp/tmpGATK -d)
echo $tmp_dir
export THEANO_FLAGS="base_compiledir=$tmp_dir"

# Load vars
script=$1
params=$2
dataset=$3
include_dir=$4
evaluate_par=$5
log=$6


# Load modules
module load apps/java-17
module load apps/R-4.1.2
module load apps/singularity

# Go to benchmark folder and run
cd /myPath/CNVbenchmarkeR2/
Rscript $script $params $dataset $include_dir $evaluate_par > $log 2>&1

# Remove tmp dr
rm -r $tmp_dir
