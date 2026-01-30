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


#Full blast results
full_blast_res <- read_tsv(paste0(opt$input,"_blockparse_targets_blast_out_btop.tsv"),
                                     col_names = c("qseqid",
                                                   "sseqid", 
                                                   "pident", 
                                                   "nident",
                                                   "length",
                                                   "mismatch", 
                                                   "gapopen",
                                                   "qstart", 
                                                   "qend", 
                                                   "qlen", 
                                                   "sstart", 
                                                   "send", 
                                                   "positive", 
                                                   "evalue",
                                                   "bitscore",
                                                   "btop",
                                                   "stitle"))


#"Hits" are those with fewer than 5 mismatches
full_blast_res_query_hit_count <- full_blast_res %>% dplyr::filter(positive > 20) %>% group_by(qseqid) %>% dplyr::summarise(count=n()) %>% ungroup()
full_blast_res_first_pass <- full_blast_res_query_hit_count %>% dplyr::filter(count < 500)

#First pass pre-filter blast results
blast_res <- full_blast_res %>% dplyr::filter(qseqid %in% full_blast_res_first_pass$qseqid)

species_to_design_for <- opt$input %>% gsub(pattern="_", replacement=" ", .)

genera_to_design_for <-  stringr::word(opt$input,1,sep="_")

#define Phocaeicola species which were wrongly named as Bacteroides in the past:
Phocaeicola_species <- c("Bacteroides_barnesiae", 
                         "Bacteroides_caecicola", 
                         "Bacteroides_caecigallinarum", 
                         "Bacteroides_chinchillae", 
                         "Bacteroides_coprocola", 
                         "Bacteroides_coprophilus", 
                         "Bacteroides_dorei", 
                         "Bacteroides_gallinaceum", 
                         "Bacteroides_massiliensis", 
                         "Bacteroides_paurosaccharolyticus", 
                         "Bacteroides_plebeius", 
                         "Bacteroides_salanitronis", 
                         "Bacteroides_sartorii", 
                         "Bacteroides_vulgatus")


Phocaeicola_mask <- "Bacteroides_barnesiae|Bacteroides_caecicola|Bacteroides_caecigallinarum|Bacteroides_chinchillae|Bacteroides_coprocola|Bacteroides_coprophilus|Bacteroides_dorei|Bacteroides_gallinaceum|Bacteroides_massiliensis|Bacteroides_paurosaccharolyticus|Bacteroides_plebeius|Bacteroides_salanitronis|Bacteroides_sartorii|Bacteroides_vulgatus"
Phocaeicola_mask_fmt <- "Bacteroides barnesiae|Bacteroides caecicola|Bacteroides caecigallinarum|Bacteroides chinchillae|Bacteroides coprocola|Bacteroides coprophilus|Bacteroides dorei|Bacteroides gallinaceum|Bacteroides massiliensis|Bacteroides paurosaccharolyticus|Bacteroides plebeius|Bacteroides salanitronis|Bacteroides sartorii|Bacteroides vulgatus"


#Define Cutibacterium species wrongly classified as Propionibacterium in the past:
Cutibacterium_species <- c("Propionibacterium_acnes", 
                           "Propionibacterium_avidum", 
                           "Propionibacterium_granulosum",
                           "Propionibacterium namnetense")

Cutibacterium_mask <- "Propionibacterium_acnes|Propionibacterium_avidum|Propionibacterium_granulosum|Propionibacterium_namnetense"
Cutibacterium_mask_fmt <- "Propionibacterium acnes|Propionibacterium avidum|Propionibacterium granulosum|Propionibacterium namnetense"

species_rename_map <- c(
  "Bacteroides barnesiae" = "Phocaeicola barnesiae",
  "Bacteroides caecicola" = "Phocaeicola caecicola",
  "Bacteroides caecigallinarum" = "Phocaeicola caecigallinarum",
  "Bacteroides chinchillae" = "Phocaeicola chinchillae",
  "Bacteroides coprocola" = "Phocaeicola coprocola",
  "Bacteroides coprophilus" = "Phocaeicola coprophilus",
  "Bacteroides dorei" = "Phocaeicola dorei",
  "Bacteroides gallinaceum" = "Phocaeicola gallinaceum",
  "Bacteroides massiliensis" = "Phocaeicola massiliensis",
  "Bacteroides paurosaccharolyticus" = "Phocaeicola paurosaccharolyticus",
  "Bacteroides plebeius" = "Phocaeicola plebeius",
  "Bacteroides salanitronis" = "Phocaeicola salanitronis",
  "Bacteroides sartorii" = "Phocaeicola sartorii",
  "Bacteroides vulgatus" = "Phocaeicola vulgatus",
  "Propionibacterium acnes" = "Cutibacterium acnes",
  "Propionibacterium avidum" = "Cutibacterium avidum",
  "Propionibacterium granulosum" = "Cutibacterium granulosum",
  "Propionibacterium namnetense" = "Cutibacterium namnetense")

species_rename_map_underscore <- c(
  "Bacteroides_barnesiae" = "Phocaeicola_barnesiae",
  "Bacteroides_caecicola" = "Phocaeicola_caecicola",
  "Bacteroides_caecigallinarum" = "Phocaeicola_caecigallinarum",
  "Bacteroides_chinchillae" = "Phocaeicola_chinchillae",
  "Bacteroides_coprocola" = "Phocaeicola_coprocola",
  "Bacteroides_coprophilus" = "Phocaeicola_coprophilus",
  "Bacteroides_dorei" = "Phocaeicola_dorei",
  "Bacteroides_gallinaceum" = "Phocaeicola_gallinaceum",
  "Bacteroides_massiliensis" = "Phocaeicola_massiliensis",
  "Bacteroides_paurosaccharolyticus" = "Phocaeicola_paurosaccharolyticus",
  "Bacteroides_plebeius" = "Phocaeicola_plebeius",
  "Bacteroides_salanitronis" = "Phocaeicola_salanitronis",
  "Bacteroides_sartorii" = "Phocaeicola_sartorii",
  "Bacteroides_vulgatus" = "Phocaeicola_vulgatus",
  "Propionibacterium_acnes" = "Cutibacterium_acnes",
  "Propionibacterium_avidum" = "Cutibacterium_avidum",
  "Propionibacterium_granulosum" = "Cutibacterium_granulosum",
  "Propionibacterium_namnetense" = "Cutibacterium_namnetense")

#Use stringr to extract the species/genera in blast results
blast_res$query_species <- stringr::word(blast_res$qseqid,1,2,sep="_")

#Replace underscores
blast_res$query_species <- gsub(pattern="_", replacement=" ",blast_res$query_species)

#substitute unnecessary symbols from stitle like "[]" in "[Ruminococcus]"
blast_res$stitle <- stringi::stri_replace_all_fixed(str=blast_res$stitle, pattern=c("[","]"), 
                                                    replacement="", vectorize_all=FALSE)

#replace misnamed species
blast_res$query_species <- recode(blast_res$query_species, !!!species_rename_map)
#extract genus
blast_res$query_genera <- stringr::word(blast_res$query_species,1,sep=" ")
#replace misnamed species in stitle
blast_res$stitle <- str_replace_all(blast_res$stitle, species_rename_map)

###########
#Drop all hits that have a positive scoring match of <= 20 because this implies more than 5 mismatches between the 25 nt query to subject

blast_res_filt <- blast_res %>% dplyr::filter(positive > 20)

##Label all "on target" blast hits

blast_res_filt$on_target_blast_species <- ifelse(str_detect(string=blast_res_filt$stitle, pattern=blast_res_filt$query_species), TRUE, FALSE)
blast_res_filt$on_target_blast_genera <- ifelse(str_detect(string=blast_res_filt$stitle, pattern=blast_res_filt$query_genera), TRUE, FALSE)

blast_res_filt <- extract_species_fn(blast_res_filt)

#Split dataframe by encoding probe
bact_probe_list <- split(blast_res_filt, f=blast_res_filt$qseqid)

#################################
#Do QC for each encoding probe.
#Under the most stringent criteria, there should be 0 matches to transcripts from different bacterial species but this is usually not possible for rRNA probes.

#Prioritize probes with one partner with a blast on-target rate of >= 90% or >=99%, wchiever is possible.
##################################
#initialize vectors
bact_probes_unfiltered <- names(bact_probe_list)

num_passes_species <- vector(length=length(bact_probes_unfiltered))

num_passes_genera <- vector(length=length(bact_probes_unfiltered))

blast_hit_num_vec <- vector(length=length(bact_probes_unfiltered))

#perform QC for each probe
for (i in 1:length(bact_probes_unfiltered)){
  #number of blast hits
  blast_hit_num <- nrow(bact_probe_list[[i]])
  on_target_sum_species <- sum(bact_probe_list[[i]]$on_target_blast_species)
  on_target_sum_genera <- sum(bact_probe_list[[i]]$on_target_blast_genera)
  
  num_passes_species[i] <- on_target_sum_species
  num_passes_genera[i] <- on_target_sum_genera
  blast_hit_num_vec[i] <- blast_hit_num
}

probe_check_summary <- data.frame(probe_name=bact_probes_unfiltered, 
                                  blast_hit_num=blast_hit_num_vec,
                                  num_specific_hits_species=num_passes_species,
                                  num_specific_hits_genera=num_passes_genera)


#Short list probes based on progressively relaxed on-target rates. Ideal > 0.99

probe_check_summary$species_on_target_rate <- probe_check_summary$num_specific_hits_species / probe_check_summary$blast_hit_num

probe_check_summary$genera_on_target_rate <- probe_check_summary$num_specific_hits_genera / probe_check_summary$blast_hit_num

blast_probe_shortlists <- blast_shortlist_fn()

##First threshold of species >= 0.99 specificity (based on # of genomes hit)
blast_res_species_99_on <- blast_res_filt %>% dplyr::filter(qseqid %in% blast_probe_shortlists[["species_99_on"]])

#Identify the UNIQUE species that are represented at this threshold
blast_res_species_99_on_species_names <- blast_res_species_99_on %>% pull(query_species) %>% unique(.) 
blast_res_species_99_on_ids <- blast_res_species_99_on %>% pull(qseqid) %>% unique(.)

#Extract all the targets at species level (on or off targets from blast results)
blast_res_species_99_on_summary_df <- lapply(bact_probe_list[blast_res_species_99_on_ids], function(df){
  
  hit_species <- unique(df$subj_species)
  
  qseqid_uniq <- unique(df$qseqid)
  
  qseqid_vec <- rep(qseqid_uniq, times=length(hit_species))
  
  num_species_hits <- rep(length(hit_species),times=length(hit_species))
  
  #query molecule
  q_molecule <- sub("_[rR]ibosomal_RNA.*","", qseqid_vec)
  
  output <- data.frame(q_molecule=q_molecule,
                       qseqid=qseqid_vec,
                       subj_species=hit_species,
                       num_species_hits=num_species_hits)
  
  return(output)
  
}) %>% do.call("rbind",.)


#Second threshold of species on target rate between 90 to 98
blast_res_species_90to98_on <- blast_res_filt %>% dplyr::filter(qseqid %in% blast_probe_shortlists[["species_90to98_on"]])

#Identify the species that are represented at this threshold
blast_res_species_90to98_on_species_names <- blast_res_species_90to98_on %>% pull(query_species) %>% unique(.) 
blast_res_species_90to98_on_ids <- blast_res_species_90to98_on %>% pull(qseqid) %>% unique(.)



#Extract all the targets at species level (on or off targets from blast results)
blast_res_species_90to98_on_summary_df <- lapply(bact_probe_list[blast_res_species_90to98_on_ids], function(df){
  
  hit_species <- unique(df$subj_species)
  
  qseqid_uniq <- unique(df$qseqid)
  
  qseqid_vec <- rep(qseqid_uniq, times=length(hit_species))
  
  num_species_hits <- rep(length(hit_species),times=length(hit_species))
  
  #query molecule
  q_molecule <- sub("_[rR]ibosomal_RNA.*","", qseqid_vec)
  
  output <- data.frame(q_molecule=q_molecule,
                       qseqid=qseqid_vec,
                       subj_species=hit_species,
                       num_species_hits=num_species_hits)
  
  return(output)
  
}) %>% do.call("rbind",.)

#Third threshold of species with relaxed on target rate
#Attempt to extract probe ids for species that are still not represented in the 99_on or 90_to_98_on ids, as backup
blast_res_species_relaxed_ids <- names(bact_probe_list)[
  !names(bact_probe_list) %in% c(blast_res_species_99_on_ids,
                                 blast_res_species_90to98_on_ids)
]

blast_res_species_relaxed <- blast_res_filt %>% dplyr::filter(qseqid %in% blast_res_species_relaxed_ids)

#Extract all the targets at species level (on or off targets from blast results)
blast_res_species_relaxed_summary_df <- lapply(bact_probe_list[blast_res_species_relaxed_ids], function(df){
  
  hit_species <- unique(df$subj_species)
  
  qseqid_uniq <- unique(df$qseqid)
  
  qseqid_vec <- rep(qseqid_uniq, times=length(hit_species))
  
  num_species_hits <- rep(length(hit_species),times=length(hit_species))
  
  #query molecule
  q_molecule <- sub("_[rR]ibosomal_RNA.*","", qseqid_vec)
  
  output <- data.frame(q_molecule=q_molecule,
                       qseqid=qseqid_vec,
                       subj_species=hit_species,
                       num_species_hits=num_species_hits)
  
  return(output)
  
}) %>% do.call("rbind",.)
################################################
#For genera

blast_res_genus_99_on <- blast_res_filt %>% 
  dplyr::filter(qseqid %in% blast_probe_shortlists[["genus_99_on"]])

blast_res_genus_99_on_ids <- blast_res_genus_99_on %>% 
  pull(qseqid) %>% unique(.)

blast_res_genus_99_on_summary_df <- lapply(bact_probe_list[blast_res_genus_99_on_ids], function(df){
  
  hit_species <- unique(df$subj_species)
  
  qseqid_uniq <- unique(df$qseqid)
  
  qseqid_vec <- rep(qseqid_uniq, times=length(hit_species))
  
  num_species_hits <- rep(length(hit_species),times=length(hit_species))
  
  #query molecule
  q_molecule <- sub("_[rR]ibosomal_RNA.*","", qseqid_vec)
  
  output <- data.frame(q_molecule=q_molecule,
                       qseqid=qseqid_vec,
                       subj_species=hit_species,
                       num_species_hits=num_species_hits)
  
  return(output)
  
}) %>% do.call("rbind",.)

#select the top 10 genus specific probes based on number of species hit

if (!is.null(blast_res_genus_99_on_summary_df) && nrow(blast_res_genus_99_on_summary_df) > 0) {
  blast_res_genus_99_on_summary_df$target_genera <- stringr::word(blast_res_genus_99_on_summary_df$q_molecule,1,sep="_")
  
  blast_res_genus_99_on_summary_df_uniq <- blast_res_genus_99_on_summary_df %>% dplyr::select(-subj_species) %>% unique(.)
  
  #get top 20 probes with broadest target range (within the same genus) per q_molecule
  
  blast_res_genus_99_on_summary_df_uniq_filt <- blast_res_genus_99_on_summary_df_uniq %>% group_by(q_molecule) %>% 
    slice_max(order_by = num_species_hits, with_ties = FALSE, n = 20) %>% ungroup()
} else {
  blast_res_genus_99_on_summary_df_uniq_filt <- data.frame()
}


#################################################

#Warning: these functions should fail gracefully if input is null

######################
#99 on species probes#
######################
#the union of main probes, left_partners and right partners, unique
blast_res_species_99_on_theoretical_pairs <- label_main_and_partners(blast_res_species_99_on %>% dplyr::select(qseqid))

############################
#90 to 98 on species probes#
############################
blast_res_species_90to98_on_theoretical_pairs <- label_main_and_partners(blast_res_species_90to98_on %>% dplyr::select(qseqid))

##########################
#Relaxed species probes###
##########################
#warning: only execute if df is not null
blast_res_species_relaxed_theoretical_pairs <- label_main_and_partners(blast_res_species_relaxed %>% dplyr::select(qseqid))

####################
#99 on genus probes#
####################
#Genus probes, shortlisted for desired genera, restricted to top 20 with broadest range
#desired genera needs to be an input in the pipeline
if (!is.null(blast_res_genus_99_on_summary_df_uniq_filt) && nrow(blast_res_genus_99_on_summary_df_uniq_filt) > 0) {
  blast_res_genus_99_on_ids_fmt <- blast_res_genus_99_on_summary_df_uniq_filt  %>%
    dplyr::select(qseqid) %>% unique(.)
} else {
  blast_res_genus_99_on_ids_fmt <- data.frame(qseqid = character(0))
}

blast_res_genus_99_on_theoretical_pairs <- label_main_and_partners(blast_res_genus_99_on_ids_fmt)



# Safely extract probe IDs, handling NULL results
probe_id_list <- list(
  if(!is.null(blast_res_species_99_on_theoretical_pairs)) blast_res_species_99_on_theoretical_pairs$main_probe_id,
  if(!is.null(blast_res_species_99_on_theoretical_pairs)) blast_res_species_99_on_theoretical_pairs$partner_probe_id,
  if(!is.null(blast_res_species_90to98_on_theoretical_pairs)) blast_res_species_90to98_on_theoretical_pairs$main_probe_id,
  if(!is.null(blast_res_species_90to98_on_theoretical_pairs)) blast_res_species_90to98_on_theoretical_pairs$partner_probe_id,
  if(!is.null(blast_res_species_relaxed_theoretical_pairs)) blast_res_species_relaxed_theoretical_pairs$main_probe_id,
  if(!is.null(blast_res_species_relaxed_theoretical_pairs)) blast_res_species_relaxed_theoretical_pairs$partner_probe_id,
  if(!is.null(blast_res_genus_99_on_theoretical_pairs)) blast_res_genus_99_on_theoretical_pairs$main_probe_id,
  if(!is.null(blast_res_genus_99_on_theoretical_pairs)) blast_res_genus_99_on_theoretical_pairs$partner_probe_id
)

# Remove NULL elements and combine
probe_id_vector <- unlist(probe_id_list[!sapply(probe_id_list, is.null)])

probes_to_extract_full_blast_res <- data.frame(probe_id = probe_id_vector) %>% unique(.)

##################################
#Format full blast res refined  ##
##################################
#This processes the BLAST results from a refined set of searches in the database to extract subject hits that are shared between main and partner probe
#This is necessary because the full blast results from the entire database caps at max 500 target sequences and hence some less selective partner probes would be missing hits that they rightfully share with their main probes
#It is okay to use the less selective partner probes if the main probes are more selective.

full_blast_res_refined <- read_tsv(paste0(opt$input,"_blockparse_targets_blast_refined_out_btop.tsv"),
                                     col_names = c("qseqid",
                                                   "sseqid", 
                                                   "pident", 
                                                   "nident",
                                                   "length",
                                                   "mismatch", 
                                                   "gapopen",
                                                   "qstart", 
                                                   "qend", 
                                                   "qlen", 
                                                   "sstart", 
                                                   "send", 
                                                   "positive", 
                                                   "evalue",
                                                   "bitscore",
                                                   "btop",
                                                   "stitle"))

candidate_full_blast_res <- full_blast_res_refined %>% dplyr::filter(qseqid %in% probes_to_extract_full_blast_res$probe_id)

#Use stringr to extract the species/genera in blast results
candidate_full_blast_res$query_species <- stringr::word(candidate_full_blast_res$qseqid,1,2,sep="_")

#Replace underscores
candidate_full_blast_res$query_species <- gsub(pattern="_", replacement=" ",candidate_full_blast_res$query_species)

#substitute unnecessary symbols from stitle like "[]" in "[Ruminococcus]"
candidate_full_blast_res$stitle <- stringi::stri_replace_all_fixed(str=candidate_full_blast_res$stitle, pattern=c("[","]"), 
                                                                   replacement="", vectorize_all=FALSE)
#replace misnamed species
candidate_full_blast_res$query_species <- recode(candidate_full_blast_res$query_species, !!!species_rename_map)
#extract genus
candidate_full_blast_res$query_genera <- stringr::word(candidate_full_blast_res$query_species,1,sep=" ")
#replace misnamed species in stitle
candidate_full_blast_res$stitle <- str_replace_all(candidate_full_blast_res$stitle, species_rename_map)

#Drop all hits that have a positive scoring match of <= 20 because this implies more than 5 mismatches between the 25 nt query to subject
candidate_full_blast_res_filt <- candidate_full_blast_res %>% dplyr::filter(positive > 20)

###################################
##Label all "on target" blast hits#
###################################

candidate_full_blast_res$on_target_blast_species <- ifelse(str_detect(string=candidate_full_blast_res$stitle, pattern=candidate_full_blast_res$query_species), TRUE, FALSE)
candidate_full_blast_res$on_target_blast_genera <- ifelse(str_detect(string=candidate_full_blast_res$stitle, pattern=candidate_full_blast_res$query_genera), TRUE, FALSE)

candidate_full_blast_res <- extract_species_fn(candidate_full_blast_res)

####
candidate_full_blast_res_filt$on_target_blast_species <- ifelse(str_detect(string=candidate_full_blast_res_filt$stitle, pattern=candidate_full_blast_res_filt$query_species), TRUE, FALSE)
candidate_full_blast_res_filt$on_target_blast_genera <- ifelse(str_detect(string=candidate_full_blast_res_filt$stitle, pattern=candidate_full_blast_res_filt$query_genera), TRUE, FALSE)

candidate_full_blast_res_filt <- extract_species_fn(candidate_full_blast_res_filt)

###########################
##Summarize BLAST results #
###########################

############################
#99 on species probes#
############################

blast_res_of_partners_for_species_main_probes_99_on_summary_df <- summarize_blast_for_candidate_pairs(blast_res_species_99_on_theoretical_pairs,
                                                                                                      blast_results = candidate_full_blast_res_filt)

blast_res_of_partners_for_species_main_probes_99_on_summary_df <- extract_and_format_query_taxon(blast_res_of_partners_for_species_main_probes_99_on_summary_df)

##
blast_res_species_99_on_summary_df <- extract_and_format_query_taxon(blast_res_species_99_on_summary_df)

# For less stringent main probes (90 to 98% specificity in BLAST results)
############################
#90 to 98 on species probes#
############################

blast_res_of_partners_for_species_main_probes_90to98_on_summary_df <- summarize_blast_for_candidate_pairs(blast_res_species_90to98_on_theoretical_pairs,
                                                                                                          blast_results = candidate_full_blast_res_filt)

blast_res_of_partners_for_species_main_probes_90to98_on_summary_df <- extract_and_format_query_taxon(blast_res_of_partners_for_species_main_probes_90to98_on_summary_df)

##
blast_res_species_90to98_on_summary_df <- extract_and_format_query_taxon(blast_res_species_90to98_on_summary_df)

########################
#Relaxed species probes#
########################

blast_res_of_partners_for_species_main_probes_relaxed_summary_df <- summarize_blast_for_candidate_pairs(blast_res_species_relaxed_theoretical_pairs,
                                                                                                        blast_results = candidate_full_blast_res_filt)

blast_res_of_partners_for_species_main_probes_relaxed_summary_df <- extract_and_format_query_taxon(blast_res_of_partners_for_species_main_probes_relaxed_summary_df)

blast_res_species_relaxed_summary_df <- extract_and_format_query_taxon(blast_res_species_relaxed_summary_df)


#Obtain theoretical LHS and RHS probe pairs. Can be null. 

blast_res_species_99_on_theoretical_pairs <- LHS_RHS_label_fn(blast_res_species_99_on_theoretical_pairs)

blast_res_species_90to98_on_theoretical_pairs <- LHS_RHS_label_fn(blast_res_species_90to98_on_theoretical_pairs)

blast_res_species_relaxed_theoretical_pairs <- LHS_RHS_label_fn(blast_res_species_relaxed_theoretical_pairs)

blast_res_genus_99_on_theoretical_pairs <- LHS_RHS_label_fn(blast_res_genus_99_on_theoretical_pairs)

######################################################
#Parse blast traceback results for LHS and RHS probes#
######################################################

blast_res_species_99_on_for_candidates <- parse_LHS_and_RHS_matches(full_blast_res=candidate_full_blast_res, 
                                                                    probe_pairs = blast_res_species_99_on_theoretical_pairs)

blast_res_species_90to98_on_for_candidates <- parse_LHS_and_RHS_matches(full_blast_res=candidate_full_blast_res, 
                                                                        probe_pairs = blast_res_species_90to98_on_theoretical_pairs)

blast_res_species_relaxed_for_candidates <- parse_LHS_and_RHS_matches(full_blast_res=candidate_full_blast_res, 
                                                                      probe_pairs = blast_res_species_relaxed_theoretical_pairs)


blast_res_genus_99_on_for_candidates <- parse_LHS_and_RHS_matches(full_blast_res=candidate_full_blast_res, 
                                                                  probe_pairs = blast_res_genus_99_on_theoretical_pairs)

#################################
#Get candidate LHS and RHS pairs#
#################################
#########
blast_res_species_99_on_candidate_pairs <- find_shared_targets(blast_res_list = blast_res_species_99_on_for_candidates,
                                                               probe_pairs = blast_res_species_99_on_theoretical_pairs)


blast_res_species_90to98_on_candidate_pairs <- find_shared_targets(blast_res_list =blast_res_species_90to98_on_for_candidates,
                                                                   probe_pairs=blast_res_species_90to98_on_theoretical_pairs)

blast_res_species_relaxed_candidate_pairs <- find_shared_targets(blast_res_list = blast_res_species_relaxed_for_candidates,
                                                                 probe_pairs = blast_res_species_relaxed_theoretical_pairs)


blast_res_genus_99_on_candidate_pairs <- find_shared_targets(blast_res_list=blast_res_genus_99_on_for_candidates,
                                                             probe_pairs=blast_res_genus_99_on_theoretical_pairs)


#write out the files.

# Helper function to safely write results
safe_write <- function(data, filepath) {
  if (!is.null(data) && is.data.frame(data) && nrow(data) > 0) {
    write_tsv(data, file = filepath)
  } else {
    # Write empty file with header if data is NULL/empty
    write_tsv(data.frame(), file = filepath)
  }
}

safe_write(blast_res_species_99_on_candidate_pairs$paired_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_99_candidate_paired_hits.tsv")))

safe_write(blast_res_species_99_on_candidate_pairs$paired_genome_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_99_candidate_paired_genome_hits.tsv")))

safe_write(blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary_fmt, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_99_candidate_hit_summary.tsv")))

####
safe_write(blast_res_species_90to98_on_candidate_pairs$paired_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_90to98_candidate_paired_hits.tsv")))

safe_write(blast_res_species_90to98_on_candidate_pairs$paired_genome_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_90to98_candidate_paired_genome_hits.tsv")))

safe_write(blast_res_species_90to98_on_candidate_pairs$probe_pair_hit_summary_fmt, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_90to98_candidate_hit_summary.tsv")))

####
safe_write(blast_res_species_relaxed_candidate_pairs$paired_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_relaxed_candidate_paired_hits.tsv")))

safe_write(blast_res_species_relaxed_candidate_pairs$paired_genome_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_relaxed_candidate_paired_genome_hits.tsv")))

safe_write(blast_res_species_relaxed_candidate_pairs$probe_pair_hit_summary_fmt, 
          file.path(output_dir,paste0(opt$input,"_blast_res_species_relaxed_candidate_hit_summary.tsv")))

###
safe_write(blast_res_genus_99_on_candidate_pairs$paired_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_genus_99_candidate_paired_hits.tsv")))

safe_write(blast_res_genus_99_on_candidate_pairs$paired_genome_hits_by_probe, 
          file.path(output_dir,paste0(opt$input,"_blast_res_genus_99_candidate_paired_genome_hits.tsv")))


safe_write(blast_res_genus_99_on_candidate_pairs$probe_pair_hit_summary_fmt, 
          file.path(output_dir,paste0(opt$input,"_blast_res_genus_99_candidate_hit_summary.tsv")))


#write out another file that just has the probe ids to extract the fasta sequences for thermoanalysis
# Safely extract probe IDs, handling NULL results
probe_id_list <- list(
  if(!is.null(blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary)) blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary$main_probe_id,
  if(!is.null(blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary)) blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary$partner_probe_id,
  if(!is.null(blast_res_species_90to98_on_candidate_pairs$probe_pair_hit_summary)) blast_res_species_90to98_on_candidate_pairs$probe_pair_hit_summary$main_probe_id,
  if(!is.null(blast_res_species_90to98_on_candidate_pairs$probe_pair_hit_summary)) blast_res_species_90to98_on_candidate_pairs$probe_pair_hit_summary$partner_probe_id,
  if(!is.null(blast_res_species_relaxed_candidate_pairs$probe_pair_hit_summary)) blast_res_species_relaxed_candidate_pairs$probe_pair_hit_summary$main_probe_id,
  if(!is.null(blast_res_species_relaxed_candidate_pairs$probe_pair_hit_summary)) blast_res_species_relaxed_candidate_pairs$probe_pair_hit_summary$partner_probe_id,
  if(!is.null(blast_res_genus_99_on_candidate_pairs$probe_pair_hit_summary)) blast_res_genus_99_on_candidate_pairs$probe_pair_hit_summary$main_probe_id,
  if(!is.null(blast_res_genus_99_on_candidate_pairs$probe_pair_hit_summary)) blast_res_genus_99_on_candidate_pairs$probe_pair_hit_summary$partner_probe_id
)

# Remove NULL elements and combine
probe_id_vector <- unlist(probe_id_list[!sapply(probe_id_list, is.null)])

all_candidate_probe_ids <- if(length(probe_id_vector) > 0) {
  data.frame(probe_id = probe_id_vector) %>% unique(.)
} else {
  data.frame(probe_id = character(0))
}

#This is an important file. It will feed into bbmap filterbyname.sh and then thermoanalysis
write_tsv(all_candidate_probe_ids, file=file.path(output_dir,paste0(opt$input,"_all_candidate_probe_ids.tsv")),col_names = FALSE)

