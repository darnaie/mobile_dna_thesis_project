#!/bin/bash
#SBATCH --job-name=extracting_flanks
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7:30:00
#SBATCH --mem=16G
#SBATCH --output=extracting_flanks_%j.out
#SBATCH --error=extracting_flanks_%j.err
#SBATCH --mail-type=END,FAIL

module load devel/python/3.10.0

python /net/bq-storage/ag-khedkar/Sofia/project_folder/extract_flanks.py

