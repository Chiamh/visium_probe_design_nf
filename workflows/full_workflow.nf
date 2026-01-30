/*
========================================================================================
    Help messages and warnings
========================================================================================
*/

def helpMessage() {
	// adapted from nf-core
	
	log.info"""
	Usage for main workflow:
	The typical command for running the pipeline is as follows:
	nextflow run main.nf
		--taxa_list PATH TO LIST OF SPECIES TO DESIGN PROBES FOR
		--genome_metadata PATH TO NCBI METADATA FOR SPECIES REPRESENTATIVE GENOMES
		--blastdb PATH TO BLAST DB FOLDER
		
	
	Input arguments:
		--taxa_list				Path to an input text file with each line containing a single microbial species, with spaces separated by underscores e.g. "Cutibacterium_acnes"
		--genome_metadata		Path to a tab separated file of NCBI metadata for species representative genome assemblies
		--blastdb				Path to the blast database to check for probe specificity
	Output arguments:
		--outdir				The output directory where the results will be saved [Default: ./pipeline_results]
		--tracedir				The directory where nextflow logs will be saved [Default: ./pipeline_results/pipeline_info]
	Others:
		--help					Display this help message
	"""
}

if (params.help){
	helpMessage()
	exit 0
}


if (!params.taxa_list){
	helpMessage()
	log.info"""
	[Error] Please specify input species for probe design
	""".stripIndent()
	exit 0
}

if (!params.genome_metadata){
	helpMessage()
	log.info"""
	[Error] Please include NCBI reference genome metadata. File format is described in README.txt
	""".stripIndent()
	exit 0
}




/*
========================================================================================
    Define channels for input taxa
========================================================================================
*/

//emits ['Staphylococcus', 'Staphylococcus_aureus'] for example
Channel
	.fromPath(params.taxa_list)
	.splitText()
	.map { line ->
		def species = line.trim()
		def genus = species.split('_', 2)[0]
		tuple(genus, species)
		}.set { ch_taxa }
		
/*
==========================================================================================
    Include modules
==========================================================================================
*/
include { DOWNLOAD_TARGETS } from '../modules/download_target_seqs.nf'
include { EXTRACT_SEQS } from '../modules/extract_target_seqs.nf'
include { PROBE_MINING } from '../modules/probemining.nf'
include { BLAST } from '../modules/blast.nf'
include { PROBE_FILTERING } from '../modules/probe_filtering.nf'
include { BLAST_REFINED } from '../modules/blast_refined.nf'
include { PROBE_FILTERING_REFINED } from '../modules/probe_filtering_refined.nf'
include { THERMOANALYSIS } from '../modules/thermoanalysis.nf'
include { PROBE_SHORTLIST } from '../modules/probe_shortlist.nf'
/*
========================================================================================
    Workflow
========================================================================================
*/

workflow FULL {	
	// 1. Download representative genomes and extract target sequences to design probes against
	
	DOWNLOAD_TARGETS( params.genome_metadata, ch_taxa )
	EXTRACT_SEQS( DOWNLOAD_TARGETS.out.refsout )
	
	// 2. Design 25 nt probes for Visium assays
	PROBE_MINING( EXTRACT_SEQS.out.target_seqs )
	
	// 3. BLAST probes against human transcriptome, tRNAs, bacterial and fungal rRNAs
	BLAST( params.blastdb, PROBE_MINING.out.blockparse_out )
	
	// 4. Filter probe candidates based on blast results.
	custom_functions_ch = channel.fromPath("${projectDir}/bin/custom_probe_filtering_functions.R")
	PROBE_FILTERING( BLAST.out.blast_res, custom_functions_ch )
	
	// 4b. Blast probe candidates with refined database
	ch_blast_refined_inputs = PROBE_MINING.out.blockparse_out.join(PROBE_FILTERING.out.tsv_files)
	BLAST_REFINED( params.blastdb, ch_blast_refined_inputs )
	
	// 4c. Filter probe candidates based refined blast results
	ch_probe_filtering_refined_inputs = BLAST.out.blast_res.join(BLAST_REFINED.out.blast_res)
	PROBE_FILTERING_REFINED( ch_probe_filtering_refined_inputs, custom_functions_ch )
	
	// 5. Analyze thermal properties of probes
	// Combine probe_ids with blockparse fasta: [species, blockparse.fa, probe_ids.tsv, blast_res*.tsv]
	
	ch_probe_thermoanalysis_inputs = PROBE_MINING.out.blockparse_out
		.join(PROBE_FILTERING_REFINED.out.tsv_files)
		.map { species, blockparse_fa, probe_ids, blast_res_files ->
			tuple(species, blockparse_fa, probe_ids)
		}
	
	THERMOANALYSIS( ch_probe_thermoanalysis_inputs )
	
	// 6. Shortlist probes
	// Recombine all files: [species, blockparse.fa, primer3.tsv, blast_res*.tsv]

	ch_probe_shortlist_inputs = PROBE_MINING.out.blockparse_out
		.join(THERMOANALYSIS.out.primer3_results)
		.join(PROBE_FILTERING_REFINED.out.tsv_files.map { species, probe_ids, blast_res -> tuple(species, blast_res) })
	PROBE_SHORTLIST( ch_probe_shortlist_inputs, custom_functions_ch )
	
}
