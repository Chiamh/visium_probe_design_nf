// For downloading representative/reference genomes from NCBI to extract target sequences for probe design

process DOWNLOAD_TARGETS {
	label "process_small"
	tag "${species}"
	publishDir "${params.outdir}/ref_genomes", mode: 'copy'
	
	input:
	path genome_metadata
	tuple val(genus), val(species)
	
	output:
	tuple val(species), path("${species}_*genomic_refseq.fna"), path("${species}_*.gff"), path("${species}_*.bed"), emit: refsout
	tuple val(species), path("${species}_*download.log"), emit: logs
	
	script:
	"""	
	
	extract_metadata.sh "${species}" "${genome_metadata}"
	
	gff2bed_batch_rRNA.R -i ${species}_list_fmt.txt
	
	
	"""
}
