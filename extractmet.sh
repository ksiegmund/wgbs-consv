#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --mail-user=youremailaddress
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50GB

cd /yourdirectory/

source /yourdirectory/setup.sh

bismark_methylation_extractor --gzip --bedGraph S${SLURM_ARRAY_TASK_ID}.multiple.deduplicated.nonCG_filtered.bam 

# echo Job in $SLURM_ARRAY_TASK_ID done
