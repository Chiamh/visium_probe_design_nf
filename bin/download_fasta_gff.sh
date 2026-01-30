#!/bin/bash

#script to download fasta and gff given some URLs

set -euo pipefail

line=$1
read -r URL SPECIES <<< "$(
	awk -F'\t' '{print $13, $7}' <<< "$line"
)"

ID=$(basename "$URL")

#Create log file in working directory if it does not already exist.
[ -f "$SPECIES"_download.log ] || touch "$SPECIES"_download.log


#Download FASTA file and unzip in the current working directory. Redirect std err to log file.
echo "$(date -u) Downloading from $URL/$ID for $SPECIES" >> "$SPECIES"_download.log
wget -O - "$URL"/"$ID"_genomic.fna.gz | gunzip -c > "$SPECIES"_genomic_refseq.fna 2>>"$SPECIES"_download.log
echo "$(date -u) $SPECIES FASTA FILE downloaded" >> "$SPECIES"_download.log

#Download GFF file and unzip in the current working directory. Redirect std err to log file.
wget -O - "$URL"/"$ID"_genomic.gff.gz | gunzip -c > "$SPECIES"_genomic_refseq.gff 2>>"$SPECIES"_download.log
echo "$(date -u) $SPECIES GFF FILE downloaded" >> "$SPECIES"_download.log

sleep 7
