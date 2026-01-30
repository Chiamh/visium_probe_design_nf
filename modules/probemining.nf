// For downloading representative/reference genomes from NCBI to extract target sequences for probe design

process PROBE_MINING {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/probe_targets", mode: 'copy'
	
	input:
	tuple val(species), path(target_fasta)
	
	output:
	tuple val(species), path("${species}_merged_blockparse_rmdup.fa"), emit: blockparse_out
	
	script:
	"""	
	
	run_oligominer.sh "${species}" "${target_fasta}"
	
	
	"""
}