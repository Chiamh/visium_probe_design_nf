process PROBE_FILTERING_REFINED {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/probe_filtering_out", mode: 'copy'
	
	input:
	tuple val(species), path(full_blast_results), path(refined_blast_results)
	each path(custom_functions)
	
	output:
	tuple val(species), path("${species}_all_candidate_probe_ids.tsv"), path("${species}_blast_res_*.tsv"), emit: tsv_files

	
	script:
	"""	
	probe_filtering_refined.R -i "${species}"
	
	"""
}
