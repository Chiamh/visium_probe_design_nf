// Analysis with primer3_py

process THERMOANALYSIS {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/probe_filtering_out", mode: 'copy'
	
	input:
	tuple val(species), path(blockparse_fasta), path(probe_ids_tsv)
	
	output:
	tuple val(species), path("${species}_primer3_thermoanalysis.tsv"), emit: primer3_results
	
	script:
	"""	
	filterbyname.sh -Xmx3g in="${blockparse_fasta}" out="${species}"_blockparse_candidates.fa names="${species}"_all_candidate_probe_ids.tsv include=t fixjunk
	
	multifasta_to_primer3.py "${species}"_blockparse_candidates.fa --mv 165 --dv 0 > "${species}"_primer3_thermoanalysis.tsv
	
	"""
}
