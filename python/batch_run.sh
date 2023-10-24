#!/bin/bash
#SBATCH --job-name=xlms_datasets
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --array=0-59%10
#SBATCH --gres=gpu:v100-sxm2:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=4GB
#SBATCH --time=04:00:00
#SBATCH --output=out/array_%A_%a.out
#SBATCH --error=err/array_%A_%a.err


# srun ./
# echo $SLURM_ARRAY_TASK_ID
. ~/.bash_profile

((i = SLURM_ARRAY_TASK_ID % 10))

./run_dataset_i.sh $i unweighted_pdf_mode



