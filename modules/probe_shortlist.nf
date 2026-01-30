process PROBE_SHORTLIST {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/probe_shortlist_out", mode: 'copy'
	
	input:
	tuple val(species), path(blockparse_fasta), path(primer3_res), path(tsv_files)
	each path(custom_functions)
	
	output:
	tuple val(species), path("${species}_{LHS,RHS}_probes_*.tsv"), path("${species}_combined_probe_candidate_summary.tsv"), emit: results

	
	script:
	"""	
	probe_shortlisting.R -i "${species}"
	
	"""
}
