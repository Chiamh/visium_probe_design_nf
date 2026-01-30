// For downloading representative/reference genomes from NCBI to extract target sequences for probe design

process EXTRACT_SEQS {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/ref_genomes", mode: 'copy'
	
	input:
	tuple val(species), path(ref_genome), path(species_gff), path(species_bed)
	
	output:
	tuple val(species), path("${species}_targets_rmdup.fa"), emit: target_seqs
	
	script:
	"""	
	
	fasta_from_bed_v2.sh "${species}"
	
	
	"""
}