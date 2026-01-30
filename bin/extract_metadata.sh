#!/bin/bash
set -euo pipefail

species=${1}
ncbi_metadata=${2}

if grep -F "${species}" "${ncbi_metadata}" > ${species}_ref_genome_metadata.tsv; then
	
	cut -f7 ${species}_ref_genome_metadata.tsv > ${species}_list_fmt.txt
	
	while IFS=$'\t' read -r line; do
		bash download_fasta_gff.sh "$line"
	done < ${species}_ref_genome_metadata.tsv
else
	echo "No genomes found for ${species}" >> "${species}_download.log"
fi




