#!/bin/bash
#SBATCH --job-name=C11_TE_lib_general
#SBATCH --output=C11_TE_lib_general_output_%j
#SBATCH --error=C11_TE_lib_general_error_%j
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=100:00:00

INDIR=/home/igomez/genome_assembly_helodes_thesis/C0_clean_assemblies
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library

    #Load conda source and earlgrey package
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate earlgrey
    #Load genom
    genome_pth=$INDIR/4JMC19C.fasta

    #mkdir -p $OUTDIR/4JMC19C
    cd $OUTDIR/4JMC19C

    #Run earlgrey
    earlGrey -g $genome_pth -s 4JMC19C_1 -l /home/igomez/genome_assembly_helodes_thesis/C1_TE_library/TE_database_personalized/plant_manual_rep_libary.fasta -o $OUTDIR/4JMC19C -t 32
    conda deactivate

