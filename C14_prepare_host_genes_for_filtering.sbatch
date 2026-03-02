#!/bin/bash
#SBATCH --job-name=C14_prepare_blast
#SBATCH --output=C14_prepare_blast_%j.out
#SBATCH --error=C14_prepare_blast_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=100:00:00


### Check if sequences from the EarlGrey library are in fact host genes
### By blasting the library against the nr Arthropoda database
OUTDIR=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes
LIB=/home/igomez/genome_assembly_helodes_thesis/C1_TE_library/host_genes/C_helodes_lib_clstrd.fasta ### LIB from EarlGrey output
FINAL_OUTPUT="$OUTDIR/Poales_nr.fasta"

###Actiavte conda
source ~/anaconda3/etc/profile.d/conda.sh
conda activate blast_env 

# Inicializar archivo de salida final vacío
> $FINAL_OUTPUT

# Procesar cada archivo comprimido de NR

cd $OUTDIR/NR
for f in nr.*.tar.gz; do
    
    # Obtener el nombre del archivo .psq resultante (asume formato estándar)
    DB_NAME=$(basename $f .tar.gz)
   
    # Extraer la parte actual
    mkdir -p $OUTDIR/NR/$DB_NAME
    cd $OUTDIR/NR/$DB_NAME
    tar -xvzf $f

    # Aplicar blastdbcmd para extraer secuencias del taxid especificado
    cd $OUTDIR
    blastdbcmd -db NR/$DB_NAME -dbtype prot -target_only -taxids 38820 -out temp.fasta
    
    # Concatenar las secuencias al archivo final
    cat temp.fasta >> $FINAL_OUTPUT

    echo "$DB_NAME added"
    
    # Eliminar archivos temporales
    rm -f temp.fasta
    rm -r NR/$DB_NAME
done

echo "fasta file complete"

conda deactivate 

### Transform fasta to diamond
conda activate diamond
diamond makedb -p 48 --in $FINAL_OUTPUT -d Viridiplantae_diamond.db
diamond blastx -p 48 -d Viridiplantae_diamond.db -q $LIB -o C_helodes_lib_blast.tsv

### Filter the blast hits so that we keep only one hit per repeat
awk -F'\t' '!seen[$1]++' C_helodes_lib_blast.tsv > C_helodes_lib_blast_uniq.tsv
### Get the list of GenBank IDs corresponding to the hits
cat C_helodes_lib_blast_uniq.tsv | cut -f 2 > C_helodes_lib_blast_uniq_GenBankID.tsv



