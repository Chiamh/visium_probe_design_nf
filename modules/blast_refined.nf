process BLAST_REFINED {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/blast_out", mode: 'copy'
	
	input:
	path blastdb
	tuple val(species), path(blockparse_fasta), path(blast_res_first_filt)
	
	output:
	tuple val(species), path("${species}_blockparse_targets_blast_refined_out_btop.tsv"), emit: blast_res
	
	script:
	"""	
	
	awk -F '\t' 'NR>1 {print \$2}' "${blast_res_first_filt}" > "${species}"_target_list.txt
	
	blastn -num_threads $task.cpus -db "${blastdb}"/*.fasta \
	-seqidlist "${species}"_target_list.txt \
	-query "${blockparse_fasta}" \
	-task blastn-short -word_size 7 -perc_identity 84 -max_target_seqs 500 \
	-outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend qlen sstart send positive evalue bitscore btop stitle" \
	-out "${species}"_blockparse_targets_blast_refined_out_btop.tsv
	
	"""
}
