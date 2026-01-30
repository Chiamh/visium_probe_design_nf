process PROBE_FILTERING {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/probe_filtering_out", mode: 'copy'
	
	input:
	tuple val(species), path(full_blast_results)
	each path(custom_functions)
	
	output:
	tuple val(species), path("${species}_blast_res_first_pass_filter.tsv"), emit: tsv_files

	
	script:
	"""	
	probe_filtering.R -i "${species}"
	
	"""
}
