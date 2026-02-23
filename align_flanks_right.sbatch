#!/bin/bash
#SBATCH --job-name=align_flanks_right
#SBATCH --partition=single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist='cln132'
#SBATCH --cpus-per-task=20
#SBATCH --time=22:30:00
#SBATCH --mem=16G
#SBATCH --error=align_flanks_right_%j.err
#SBATCH --output=align_flanks_right_%j.out
#SBATCH --mail-type=END,FAIL

module load devel/python/3.10.0

/home/bq_slipaeva/.local/bin/vclust align \
-i /net/bq-storage/ag-khedkar/Sofia/project_folder/right_flanks.fna \
-o /net/bq-storage/ag-khedkar/Sofia/project_folder/flanks_alignment/right_ani.tsv \
--outfmt standard \
--out-ani 0.98 \
--out-qcov 1.0 \
--out-rcov 1.0