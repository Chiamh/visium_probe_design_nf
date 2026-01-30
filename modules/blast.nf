// For downloading representative/reference genomes from NCBI to extract target sequences for probe design

process BLAST {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/blast_out", mode: 'copy'
	
	input:
	path blastdb
	tuple val(species), path(blockparse_fasta)
	
	output:
	tuple val(species), path("${species}_blockparse_targets_blast_out_btop.tsv"), emit: blast_res
	
	script:
	"""	
	
	blastn -num_threads $task.cpus -db "${blastdb}"/*.fasta -query "${species}"_merged_blockparse_rmdup.fa \
	-task blastn-short -word_size 7 -perc_identity 84 -max_target_seqs 500 \
	-outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend qlen sstart send positive evalue bitscore btop stitle" \
	-out "${species}"_blockparse_targets_blast_out_btop.tsv
	
	"""
}
