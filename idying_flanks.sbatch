#!/bin/bash
#SBATCH --job-name=idying_flanks
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist='cln132'
#SBATCH --cpus-per-task=20
#SBATCH --time=22:30:00
#SBATCH --mem=16G
#SBATCH --error=idying_flanks_%j.err
#SBATCH --mail-type=END,FAIL

module load devel/python/3.10.0

python /net/bq-storage/ag-khedkar/Sofia/project_folder/idying_flanks.py
