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

safe_write(blast_res_filt, 
          file.path(output_dir,paste0(opt$input,"_blast_res_first_pass_filter.tsv")))


