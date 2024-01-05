#!/bin/bash
#SBATCH --job-name=xlms
#SBATCH --nodes=1
#SBATCH --array=0-59
#SBATCH --cpus-per-task=10
#SBATCH --mem=4GB
#SBATCH --time=04:00:00
#SBATCH --output=out/array_%A_%a.out
#SBATCH --error=err/array_%A_%a.err


# srun ./
# echo $SLURM_ARRAY_TASK_ID
. ~/.bash_profile

((i = SLURM_ARRAY_TASK_ID % 10))
((p = SLURM_ARRAY_TASK_ID / 10))
mu_strategy=$1

echo $mu_strategy

./run_1S_dataset_i.sh $i unweighted_pdf_mode $p $mu_strategy
# ./run_dataset_i.sh $i unweighted_pdf_mode $p $mu_strategy


