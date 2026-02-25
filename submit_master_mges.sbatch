#!/bin/bash
#SBATCH --job-name=master_mges
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=05:30:00
#SBATCH --mem=32G
#SBATCH --output=master_mges_%j.out
#SBATCH --error=master_mges_%j.err
#SBATCH --mail-type=END,FAIL

module load python/3.9

python /net/bq-storage/ag-khedkar/Sofia/project_folder/gff3_and_ffn_file_processing.py
