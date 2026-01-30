#!/usr/bin/env Rscript

# GFF to BED Converter - Generalizable RNA Feature Extraction
# This script converts GFF files to BED format for specified RNA features
# Date: 2026-01-18

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(biomartr)
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input species list file (TSV with 'species' column)", metavar = "FILE"),
  make_option(c("-d", "--gff_dir"), type = "character", default = ".",
              help = "Directory containing GFF files [default: %default]", metavar = "DIR"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Output directory for BED files [default: same as gff_dir]", metavar = "DIR"),
  make_option(c("-s", "--gff_suffix"), type = "character", default = "_genomic_refseq.gff",
              help = "Suffix for GFF filenames [default: %default]", metavar = "SUFFIX"),
  make_option(c("-b", "--bed_suffix"), type = "character", default = "_targets.bed",
              help = "Suffix for output BED filenames [default: %default]", metavar = "SUFFIX"),
  make_option(c("-f", "--feature_type"), type = "character", default = "rRNA",
              help = "Feature type to extract from GFF (e.g., rRNA, tRNA, CDS) [default: %default]", metavar = "TYPE"),
  make_option(c("-p", "--pattern"), type = "character", default = "16S|23S|18S|28S",
              help = "Pattern to filter features by (regex) [default: %default]", metavar = "PATTERN"),
  make_option(c("-a", "--attribute_field"), type = "character", default = "product",
              help = "GFF attribute field to extract for BED name [default: %default]", metavar = "FIELD"),
  make_option(c("--warn_missing"), action = "store_true", default = TRUE,
              help = "Warn if no features detected in GFF file [default: %default]"),
  make_option(c("--quiet"), action = "store_true", default = FALSE,
              help = "Suppress progress messages [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input species list file must be specified with -i/--input", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call. = FALSE)
}

# Set output directory to gff_dir if not specified
if (is.null(opt$output_dir)) {
  opt$output_dir <- opt$gff_dir
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  if (!opt$quiet) message("Created output directory: ", opt$output_dir)
}

# Read input species list
input_species <- read_tsv(opt$input, col_names = "species", show_col_types = FALSE)

if (!opt$quiet) {
  message(sprintf("Processing %d species from %s", nrow(input_species), opt$input))
  message(sprintf("Feature type: %s, Filter pattern: %s", opt$feature_type, opt$pattern))
}

# Generalizable function to make BED files
make_bed_files <- function(species, gff_dir, output_dir, gff_suffix, bed_suffix, 
                           feature_type, pattern, attribute_field, warn_missing, quiet) {
  
  gff_file <- file.path(gff_dir, paste0(species, gff_suffix))
  
  # Check if GFF file exists
  if (!file.exists(gff_file)) {
    warning(sprintf("GFF file not found for %s: %s", species, gff_file))
    return(invisible(NULL))
  }
  
  # Read GFF file
  tryCatch({
    gff_df <- read_gff(file = gff_file)
  }, error = function(e) {
    warning(sprintf("Error reading GFF file for %s: %s", species, e$message))
    return(invisible(NULL))
  })
  
  # Filter by feature type
  filtered_gff <- gff_df %>% dplyr::filter(type == feature_type)
  
  # Apply pattern filter if specified
  if (!is.null(pattern) && pattern != "") {
    filtered_gff <- filtered_gff %>% 
      dplyr::filter(str_detect(.$attribute, pattern))
  }
  
  # Check if any features were found
  if (nrow(filtered_gff) == 0) {
    if (warn_missing) {
      warning(sprintf("No %s features matching pattern '%s' found for %s", 
                      feature_type, pattern, species))
    }
    return(invisible(NULL))
  }
  
  # Extract attribute field for BED name
  regex_pattern <- paste0(attribute_field, "\\=(.*?)(;|$)")
  attr_match <- stringr::str_match(filtered_gff$attribute, regex_pattern)[, 2]
  
  # Handle cases where attribute field is not found
  filtered_gff$mod_attr <- ifelse(is.na(attr_match), 
                                  paste0(species, "_", feature_type, "_", seq_len(nrow(filtered_gff))),
                                  attr_match)
  
  # Clean up attribute names
  filtered_gff$mod_attr <- gsub(pattern = " ", replacement = "_", x = filtered_gff$mod_attr)
  filtered_gff$mod_attr <- paste0(species, "_", filtered_gff$mod_attr)
  
  # Create BED format output
  output_bed <- filtered_gff %>% 
    dplyr::select(c(seqid, start, end, mod_attr, score, strand))
  
  # Convert from 1-based coords (GFF) to 0-based coords (BED)
  # Reference: https://tidyomics.com/blog/2018/12/09/2018-12-09-the-devil-0-and-1-coordinate-system-in-genomics/
  output_bed$start <- output_bed$start - 1
  
  # Write BED file
  bed_file <- file.path(output_dir, paste0(species, bed_suffix))
  write_tsv(output_bed, file = bed_file, col_names = FALSE)
  
  if (!quiet) {
    message(sprintf("Created BED file for %s: %d features written to %s", 
                    species, nrow(output_bed), bed_file))
  }
  
  return(invisible(output_bed))
}

# Process all species
if (!opt$quiet) message("\nProcessing species...")

success_count <- 0
error_count <- 0

for (i in input_species$species) {
  if (!opt$quiet) message(sprintf("Processing: %s", i))
  
  result <- tryCatch({
    make_bed_files(
      species = i,
      gff_dir = opt$gff_dir,
      output_dir = opt$output_dir,
      gff_suffix = opt$gff_suffix,
      bed_suffix = opt$bed_suffix,
      feature_type = opt$feature_type,
      pattern = opt$pattern,
      attribute_field = opt$attribute_field,
      warn_missing = opt$warn_missing,
      quiet = opt$quiet
    )
    success_count <<- success_count + 1
    TRUE
  }, error = function(e) {
    error_count <<- error_count + 1
    warning(sprintf("Failed to process %s: %s", i, e$message))
    FALSE
  })
}

# Summary
if (!opt$quiet) {
  message("\n=== Processing Complete ===")
  message(sprintf("Total species: %d", nrow(input_species)))
  message(sprintf("Successfully processed: %d", success_count))
  message(sprintf("Errors: %d", error_count))
}
