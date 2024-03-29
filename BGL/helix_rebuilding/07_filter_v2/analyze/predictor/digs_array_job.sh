#!/bin/bash
#SBATCH -J predict_error
#SBATCH -p backfill
#SBATCH --mem=8g
#SBATCH -o digs.log
#SBATCH -e digs.log

echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date)"
CMD=$(head -n $SLURM_ARRAY_TASK_ID tasks | tail -1)
exec ${CMD}

#####command:
#####sbatch -a 1-$(cat tasks|wc -l) digs_array_job.sh
