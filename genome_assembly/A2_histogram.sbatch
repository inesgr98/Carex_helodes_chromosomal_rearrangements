#!/usr/bin/bash
#SBATCH --job-name=A2_generate_histogram
#SBATCH --output=A2_generate_histogram_output_%j
#SBATCH --error=A2_generate_histogram_error_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=standard
#SBATCH --time=2:00:00


source ~/anaconda3/etc/profile.d/conda.sh
conda activate purge_duplicates
cd ~/genome_assembly_helodes_thesis/A2_purge_dups/1JMC18

    hist_plot.py -c cutoffs PB.stat cov_unpurged_hist.png
