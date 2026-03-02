#!/bin/bash
#SBATCH --job-name=C12_TE_lib_conc
#SBATCH --output=C12_TE_lib_conc_output_%j
#SBATCH --error=C12_TE_lib_conc_error_%j
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=100:00:00

INDIR=/home/igomez/genome_assembly_helodes_thesis/C0_clean_assemblies
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library

    #Load conda source and earlgrey package
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate earlgrey
    #Load genom
    genome_pth=$INDIR/1JMC18.fasta

    mkdir -p $OUTDIR/1JMC18
    cd $OUTDIR/1JMC18

    #Run earlgrey
    earlGrey -g $genome_pth -s 1JMC18r1_conc -l /home/igomez/genome_assembly_helodes_thesis/C1_TE_library/4JMC19C/4JMC19C_r1_EarlGrey/4JMC19C_r1_Curated_Library/4JMC19C_r1_combined_library.fasta -o $OUTDIR/1JMC18 -t 32
    conda deactivate

