#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --mail-user=youremailaddress
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=main
#SBATCH --mem=150GB

cd /yourdirectory/

# Specifiy where your setup file is (help to locate the Bismark program)
source /yourdirectory/setup.sh

deduplicate_bismark --bam -p --multiple S${SLURM_ARRAY_TASK_ID}.bam S${SLURM_ARRAY_TASK_ID}a.bam S${SLURM_ARRAY_TASK_ID}b.bam

filter_non_conversion -p  S${SLURM_ARRAY_TASK_ID}.multiple.deduplicated.bam --threshold 1  

# echo Job in $SLURM_ARRAY_TASK_ID done
