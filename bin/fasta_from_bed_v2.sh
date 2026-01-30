#!/bin/bash

#This script is downstream of gff2bed_batch.R
#It is used to extract target sequences defined in a bed file (like rRNAs and hsekeeping genes) from the *.genomic_refseq.fna
#requires bedtools and seqkit from the bioinfo conda environment or equivalent

#The input is a single species name as a string, e.g. Bacteroides_ovatus
#For convenience, run this in the same wd as the fna and bed files.
#Example usage: bash fasta_from_bed_v2.sh Staphylococcus_aureus

set -euo pipefail

# Check if species name is provided
if [ -z "$1" ]; then
    echo "Error: No species name provided"
    echo "Usage: bash fasta_from_bed_v2.sh <species_name>"
    echo "Example: bash fasta_from_bed_v2.sh Staphylococcus_aureus"
    exit 1
fi

SPECIES=$1

echo "extracting target sequences for $SPECIES"
bedtools getfasta -fi "$SPECIES"_*genomic_refseq.fna -bed "$SPECIES"_*targets.bed -s -name > "$SPECIES"_temp.fa

picard NormalizeFasta I="$SPECIES"_temp.fa O="$SPECIES"_targets.fa LINE_LENGTH=70

rm "$SPECIES"_temp.fa
echo "extraction done for $SPECIES"

#Remove duplicate target sequences
seqkit rmdup -s "$SPECIES"_targets.fa > "$SPECIES"_targets_rmdup.fa
#Remove problematic characters from fasta headers
sed -i 's/(+)//g; s/(-)//g; s/(/_/g; s/)/_/g' "$SPECIES"_targets_rmdup.fa
