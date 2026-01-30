#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(optparse)
})

source("custom_probe_filtering_functions.R")

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Species name separated with underscore e.g. Cutibacterium_acnes"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Output directory", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  stop("ERROR: --input is required (e.g. Cutibacterium_acnes)")
}

# Set output directory
output_dir <- opt$output_dir

#thermoanalysis results
probe_primer3 <- read_tsv(paste0(opt$input,"_primer3_thermoanalysis.tsv"))

probe_primer3$target_self_Tm_difference <- probe_primer3$MeltTemp - probe_primer3$Hairpin_Tm_C

probe_primer3_passed <- probe_primer3 %>% dplyr::filter(Hairpin_Tm_C <= 45 & Homodimer_Tm_C <= 45 & target_self_Tm_difference >= 10)


#filter probes (blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary_fmt) for only those that have passed thermoanalysis
blast_res_species_99_on_candidate_hit_summary_passed <- read_tsv(paste0(opt$input,"_blast_res_species_99_candidate_hit_summary.tsv")) %>% filter_by_thermoanalysis(.)
if (!is.null(blast_res_species_99_on_candidate_hit_summary_passed)) {
  blast_res_species_99_on_candidate_hit_summary_passed$category <- "species_candidate"
}

blast_res_species_90to98_on_candidate_hit_summary_passed <- read_tsv(paste0(opt$input,"_blast_res_species_90to98_candidate_hit_summary.tsv")) %>% filter_by_thermoanalysis(.)
if (!is.null(blast_res_species_90to98_on_candidate_hit_summary_passed)) {
  blast_res_species_90to98_on_candidate_hit_summary_passed$category <- "species_candidate"
}

blast_res_species_relaxed_candidate_hit_summary_passed <- read_tsv(paste0(opt$input,"_blast_res_species_relaxed_candidate_hit_summary.tsv")) %>% filter_by_thermoanalysis(.)
if (!is.null(blast_res_species_relaxed_candidate_hit_summary_passed)) {
  blast_res_species_relaxed_candidate_hit_summary_passed$category <- "species_candidate"
}

blast_res_genus_99_on_candidate_hit_summary_passed <- read_tsv(paste0(opt$input,"_blast_res_genus_99_candidate_hit_summary.tsv")) %>% filter_by_thermoanalysis(.)
if (!is.null(blast_res_genus_99_on_candidate_hit_summary_passed)) {
  blast_res_genus_99_on_candidate_hit_summary_passed$category <- "genus_candidate"
}

#combine the blast res summaries with sequence info
# Only bind non-NULL dataframes
candidate_list <- list(
  blast_res_species_99_on_candidate_hit_summary_passed,
  blast_res_species_90to98_on_candidate_hit_summary_passed,
  blast_res_species_relaxed_candidate_hit_summary_passed,
  blast_res_genus_99_on_candidate_hit_summary_passed
)
candidate_list <- candidate_list[!sapply(candidate_list, is.null)]

if (length(candidate_list) > 0) {
  combined_probe_candidate_summary <- dplyr::bind_rows(candidate_list) %>% unique(.)
  combined_probe_candidate_out <- append_Visium_adapters(combined_probe_candidate_summary)
} else {
  message("Warning: No probes passed all filters")
  combined_probe_candidate_summary <- NULL
  combined_probe_candidate_out <- NULL
}

# Helper function to safely write results
safe_write <- function(data, filepath) {
  if (!is.null(data) && is.data.frame(data) && nrow(data) > 0) {
    write_tsv(data, file = filepath)
  } else {
    # Write empty file with header if data is NULL/empty
    write_tsv(data.frame(), file = filepath)
  }
}

safe_write(if(!is.null(combined_probe_candidate_summary)) combined_probe_candidate_summary else NULL, 
          file.path(output_dir,paste0(opt$input,"_combined_probe_candidate_summary.tsv")))

safe_write(if(!is.null(combined_probe_candidate_out)) combined_probe_candidate_out$LHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_LHS_probes_before_shortlist.tsv")))

safe_write(if(!is.null(combined_probe_candidate_out)) combined_probe_candidate_out$RHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_RHS_probes_before_shortlist.tsv")))


#Shortlist probes from the combined dataframe. This requires some assumptions and should be manually checked by the user.

if (!is.null(combined_probe_candidate_out) && length(candidate_list) > 0) {

combined_probe_candidate_summary$query_species <- stringr::word(combined_probe_candidate_summary$main_probe_id,1,2,sep="_")
combined_probe_candidate_summary$query_genera <- stringr::word(combined_probe_candidate_summary$query_species,1,sep="_")

combined_probe_candidate_summary$query_molecule <- stringr::word(combined_probe_candidate_summary$main_probe_id,1,sep="::pos:")

#extract coordinates
combined_probe_candidate_summary$coordinates_main <- stringr::word(combined_probe_candidate_summary$main_probe_id,2,sep="::pos:")
combined_probe_candidate_summary$coordinates_partner <- stringr::word(combined_probe_candidate_summary$partner_probe_id,2,sep="::pos:")

combined_probe_candidate_summary$main_probe_start_coordinate <- stringr::word(combined_probe_candidate_summary$coordinates_main,1,sep="-") %>% as.numeric()
combined_probe_candidate_summary$main_probe_end_coordinate <- stringr::word(combined_probe_candidate_summary$coordinates_main,2,sep="-") %>% as.numeric()

combined_probe_candidate_summary$partner_probe_start_coordinate <- stringr::word(combined_probe_candidate_summary$coordinates_partner,1,sep="-") %>% as.numeric()
combined_probe_candidate_summary$partner_probe_end_coordinate <- stringr::word(combined_probe_candidate_summary$coordinates_partner,2,sep="-") %>% as.numeric()


#For species probes
combined_probe_candidate_summary_SPECIES <- combined_probe_candidate_summary %>% 
  dplyr::filter(category=="species_candidate")

if (nrow(combined_probe_candidate_summary_SPECIES) == 0) {
  suggested_species_probes_overlap_checked <- NULL
} else {
valid_species <- combined_probe_candidate_summary_SPECIES$query_species %>% unique(.)

suggested_species_probes_overlap_checked <- lapply(valid_species, function(species_name){
  #per "valid species"
  df_species <- combined_probe_candidate_summary_SPECIES %>% dplyr::filter(query_species == species_name)
  
  output <- select_species_probes_from_shortlist(df_species)
  
  return(output)
}) %>% do.call("rbind",.)
}


#For genus probes
combined_probe_candidate_summary_GENUS <- combined_probe_candidate_summary %>% 
  dplyr::filter(category=="genus_candidate")

if (nrow(combined_probe_candidate_summary_GENUS) == 0) {
  suggested_genus_probes_overlap_checked <- NULL
} else {
valid_genus <- combined_probe_candidate_summary_GENUS$query_genera %>% unique(.)

suggested_genus_probes_overlap_checked <- lapply(valid_genus, function(genus_name){
  #per "valid genus"
  df_genus <- combined_probe_candidate_summary_GENUS %>% dplyr::filter(query_genera == genus_name)
  
  output <- select_genus_probes_from_shortlist(df_genus)
  
  return(output)
}) %>% do.call("rbind",.)
}

#save LHS and RHS shortlisted probes

if (!is.null(suggested_species_probes_overlap_checked) && nrow(suggested_species_probes_overlap_checked) > 0) {
  shortlisted_species_probe_out <- append_Visium_adapters(suggested_species_probes_overlap_checked %>% 
                                                            dplyr::filter(shortlisted==TRUE))
} else {
  shortlisted_species_probe_out <- NULL
}

if (!is.null(suggested_genus_probes_overlap_checked) && nrow(suggested_genus_probes_overlap_checked) > 0) {
  shortlisted_genus_probe_out <- append_Visium_adapters(suggested_genus_probes_overlap_checked %>% 
                                                            dplyr::filter(shortlisted==TRUE))
} else {
  shortlisted_genus_probe_out <- NULL
}

safe_write(if(!is.null(shortlisted_species_probe_out)) shortlisted_species_probe_out$LHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_LHS_probes_shortlisted_SPECIES.tsv")))

safe_write(if(!is.null(shortlisted_species_probe_out)) shortlisted_species_probe_out$RHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_RHS_probes_shortlisted_SPECIES.tsv")))

safe_write(if(!is.null(shortlisted_genus_probe_out)) shortlisted_genus_probe_out$LHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_LHS_probes_shortlisted_GENUS.tsv")))

safe_write(if(!is.null(shortlisted_genus_probe_out)) shortlisted_genus_probe_out$RHS_all_probes else NULL, 
          file.path(output_dir,paste0(opt$input,"_RHS_probes_shortlisted_GENUS.tsv")))

} else {
  message("Warning: No combined probe candidates available for shortlisting")
  # Write empty files
  safe_write(NULL, file.path(output_dir,paste0(opt$input,"_LHS_probes_shortlisted_SPECIES.tsv")))
  safe_write(NULL, file.path(output_dir,paste0(opt$input,"_RHS_probes_shortlisted_SPECIES.tsv")))
  safe_write(NULL, file.path(output_dir,paste0(opt$input,"_LHS_probes_shortlisted_GENUS.tsv")))
  safe_write(NULL, file.path(output_dir,paste0(opt$input,"_RHS_probes_shortlisted_GENUS.tsv")))
}

