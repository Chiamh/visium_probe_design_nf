# Nextflow pipeline for Visium FFPE custom probe design for microbes

## Introduction

Chiamh/visium_probe_design_nf is a bioinformatics pipeline that takes an input list of microbial (bacterial or fungal) species and designs LHS and RHS probes for detection in Visium FFPE assays.

The pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) , a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It uses Docker containers (also compatible with Singularity) for ease of installation and computational reproducibility. 

## Pipeline summary
1. Target sequences e.g. ribosomal RNA sequences are downloaded from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) and extracted.
2. Regions from target molecules are selected using [oligominer](https://github.com/beliveau-lab/OligoMiner). These tile all eligible sequences for each target molecule.
3. Tiling oligonucleotide sequences (probes) are blasted against a database. By default, this database includes the human transcriptome, tRNAs, SSU and LSU rRNA sequences from [Silva 138.1](https://www.arb-silva.de/) 
4. Candidate probes are initially filtered to maximize species or genus specificity, and then for desired thermodynamic properties (see section on probe design parameters).
5. Non-overlapping Species- and Genus- level probes are shortlisted and the actual sequences are appended with Visium adapter sequences.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`) and add the nextflow executable to your $PATH

2. Install [`Docker`](https://docs.docker.com/engine/installation/)   

3. Clone the pipeline and refer to the help message
	```sh
	$ git clone https://github.com/Chiamh/visium_probe_design-nf
	
	$ nextflow run ./meta-omics-nf/main.nf --help
	```
* Add a custom config file which contains the paths to various pre-installed databases. Refer to the test.config file in this repo for an example. 
* Add a custom profile in the nextflow.config file, allowing you to specify the use of docker or singularity, and/or a task scheduler.  

4. Make sure all helper scripts in visium_probe_design-nf/bin have execute permissions

	```sh
	$ chmod +x ./visium_probe_design-nf/bin/*
	```

5. Run the full workflow
* Add the -bucket-dir argument if running on AWSbatch with S3 support
	```sh
	$ nextflow run ./visium_probe_design-nf/main.nf -profile docker,your_profile --taxa_list microbes_to_design.txt --outdir /path/to/results
	```
	
## Input requirements
1. Path to a text file e.g. microbes to design.txt, with each line being a microbial species of interest, with the genus and species names separated by an underscore.
<img src='/docs/input_example.png' width='100'>

## Output files

* probe_shortlist_out
	* \${SPECIES}_combined_probe_candidate_summary.tsv : Summary of BLAST results for probes for a given species. Useful for manual selection of additional probes.
	* \*_LHS_probes_before_shortlist.tsv : Left hand visium probes that were designed from target sequences for a given species. Not yet shortlisted, but included if the user wants to manually select probes.
	* \*_RHS_probes_before_shortlist.tsv : Right hand visium probes that were designed from target sequences for a given species. Not yet shortlisted.
	* \*__LHS_probes_shortlisted_SPECIES.tsv : Shortlisted species level probes (left handed). See section on probe design for how these were shortlisted.
	* \*__RHS_probes_shortlisted_SPECIES.tsv : Shortlisted species level probes (right handed). See section on probe design for how these were shortlisted.
	* \*__LHS_probes_shortlisted_GENUS.tsv : Shortlisted genus level probes (left handed). See section on probe design for how these were shortlisted.
	* \*__RHS_probes_shortlisted_GENUS.tsv : Shortlisted genus level probes (right handed). See section on probe design for how these were shortlisted.
	
	
## Probe Design Parameters
**IMPORTANT: Shortlisted probes are not always suitable! We encourage the user to select from the set of shortlisted probes, and manually check them by BLAST again.** <br>
**10X Genomics recommends selecting three probe pairs per target RNA, although this might not always be possible**
**There may also be alternative probes that can be chosen in the before_shortlisting files.** <br>
Design parameters were in part, derived from 10X Genomics's Technical Note "Custom Probe Design for Visium Spatial Gene Expression and Chromium Single Cell Gene Expression Flex" (CG000621)
1. LHS and RHS probes are 25 bp in length.
2. GC content is between 44-72% for each 25 bp probe half.
3. SNPs and mismatches at the ligation junction are avoided for "on-target" hits.
4. Matches to off-target genes have at least have five mismatches in at least one of the LHS or RHS probes to prevent efficient hybridization. 
5. An off-target with more than 75% (> 20 bp match) total complementarity to a probe is a risk for cross-hybridization, but this can be tolerated if the other partner (LHS or RHS) does not cross-hybridize to the same off-target.
6. Probe pairs do not overlap to avoid competition for the same binding sites.
7. Tm range is 48 - 77 degrees Celsius. Melting temperatures were calculated based on Visium hyb conditions which are 0% formamide and 165 mM Na+ (1X SSC buffer).
8. Self-hairpin Tm and self-dimer Tm are <= 45 degrees Celsius, which is >= 5 degrees below the hybridization temp of 50 degrees Celsius.
9. The difference between the Tm of hybridization to target and the Tm of oligo self-hairpins (secondary structure) is > 10 degrees Celsius to avoid self inhibition
10. Shortlisted species probes are based on the top 7 probe pairs which had the lowest number of off-target species by BLAST.
11. Shortlisted genus probes are based on the top 7 probe pairs which had the lowest number of off-target genera by BLAST, then sorted for the highest number of on-target species hits.

**Once again, I strongly encourage the user to check the BLAST results for the shortlisted probes in the blast_out folder, and/or manually blast them before ordering.**
	
## Contact
Minghao Chia: chia_minghao@a-star.edu.sg, chiaminghao@gmail.com
