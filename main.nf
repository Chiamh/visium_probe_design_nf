#!/usr/bin/env nextflow

/*
========================================================================================
    10X Visium FFPE Custom Probe Design
========================================================================================
    Github : 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

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


/*
========================================================================================
    Include modules
========================================================================================
*/

include { FULL } from './workflows/full_workflow.nf'

/*
========================================================================================
    Main workflow (default)
========================================================================================
*/

// this main workflow will generate all intermediate files and can be resumed with nextflow. 
workflow {
    
    FULL ()
     
}