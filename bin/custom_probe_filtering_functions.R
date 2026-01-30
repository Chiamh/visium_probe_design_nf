#########################################
#functions for probe candidate selection#
#########################################

#Extract species level ID of the subject (in the blast db). Include up to family level info because there are some uncultured bacterium
#df is something like blast_res
extract_species_fn <- function(df){
  
  df$subj_species <- stringr::word(df$stitle,-3,-1,sep=";")
  
  return(df)
  
}

### Function to get shortlists at progressive thresholds

blast_shortlist_fn <- function(df=probe_check_summary){
  
  species_99_on <-df %>% dplyr::filter(species_on_target_rate >= 0.99) %>% pull(probe_name)
  species_90to98_on <- df %>% dplyr::filter(species_on_target_rate >= 0.90 &
                                              species_on_target_rate <0.99) %>% pull(probe_name)
  
  genus_99_on <-df %>% dplyr::filter(genera_on_target_rate >= 0.99) %>% pull(probe_name)
  genus_90to98_on <-df %>% dplyr::filter(genera_on_target_rate >= 0.90 &
                                           genera_on_target_rate <0.99) %>% pull(probe_name)
  
  
  
  output <- tibble::lst(species_99_on,species_90to98_on, genus_99_on, genus_90to98_on)
  
}

#function to find all theoretical partner probe ids. Accepts a one column df with col name "qseqid"
generate_partner_probe_ids <- function(main_probe_ids){
  
  main_probe_temp <- main_probe_ids
  main_probe_temp$query_molecule <- stringr::word(main_probe_temp$qseqid,1,sep="::pos:")
  vector_of_coordinates <- stringr::word(main_probe_temp$qseqid,2,sep="::pos:")
  main_probe_temp$start_coordinate <-  stringr::word(vector_of_coordinates,1,sep="-") %>% as.numeric()
  main_probe_temp$end_coordinate <-  stringr::word(vector_of_coordinates,2,sep="-") %>% as.numeric()
  
  output_list <- list()
  for (i in 1:nrow(main_probe_temp)){
    
    main_probe_qseqid <- main_probe_temp[i, "qseqid"] %>% pull(qseqid)
    main_probe_query_mol <- main_probe_temp[i, "query_molecule"] %>% pull(query_molecule)
    main_probe_start_coord <- main_probe_temp[i, "start_coordinate"] %>% pull(start_coordinate)
    main_probe_end_coord <- main_probe_temp[i, "end_coordinate"] %>% pull(end_coordinate)
    
    left_partner_probe_query_mol <- main_probe_query_mol
    right_partner_probe_query_mol <- main_probe_query_mol
    
    left_partner_end_coord <- main_probe_start_coord -1
    left_partner_start_coord <- left_partner_end_coord - 24
    
    right_partner_start_coord <- main_probe_end_coord + 1
    right_partner_end_coord <- right_partner_start_coord + 24
    
    left_partner_qseqid <- paste0(left_partner_probe_query_mol,"::pos:",left_partner_start_coord,"-",left_partner_end_coord)
    right_partner_qseqid <- paste0(right_partner_probe_query_mol,"::pos:",right_partner_start_coord,"-",right_partner_end_coord)
    
    output_list[[i]] <- data.frame(qseqid=c(main_probe_qseqid,left_partner_qseqid,right_partner_qseqid))
    
  }
  
  output_df <- do.call("rbind", output_list)
  
  return(output_df)
  
}

#Helper function to flag out the main species probes and to display their partner probes in another column (if they exist)
#input_df consists of two columns: qseqid and probe_type, which is a character vector with elements being either: "main_probe" or "partner_probe" 
generate_probes_w_partners <- function(input_df){
  
  output <- input_df %>% as.data.frame()
  output$query_molecule <- stringr::word(output$qseqid,1,sep="::pos:")
  
  vector_of_coordinates <- stringr::word(output$qseqid,2,sep="::pos:")
  
  output$start_coordinate <-  stringr::word(vector_of_coordinates,1,sep="-") %>% as.numeric()
  output$end_coordinate <-  stringr::word(vector_of_coordinates,2,sep="-") %>% as.numeric()
  
  #For each main probe (chosen as "main" for specificity to a target), look for existence of a valid partner probe. 
  #Note: A partner probe can also be a main probe. Main probe just refers to the probes that were in the initial shortlist
  
  #Temporary dfs
  main_probe_temp <- output %>% dplyr::filter(probe_type == "main_probe")
  partner_probe_temp <- output
  
  output_list <- list()
  
  for (i in 1:nrow(main_probe_temp)){
    
    main_probe_qseqid <- main_probe_temp[i, "qseqid"]
    main_probe_query_mol <- main_probe_temp[i, "query_molecule"]
    main_probe_start_coord <- main_probe_temp[i, "start_coordinate"]
    main_probe_end_coord <- main_probe_temp[i, "end_coordinate"]
    
    #partner_probe_selection is still a tibble. Will be an empty tibble if no hits
    partner_probe_selection <- partner_probe_temp %>% 
      dplyr::filter(query_molecule == main_probe_query_mol & 
                      (start_coordinate == (main_probe_end_coord + 1) |
                         end_coordinate == (main_probe_start_coord -1)))
    
    if(nrow(partner_probe_selection>=1)){
      partner_probe_selection$main_probe_id <- main_probe_qseqid
      partner_probe_selection$main_probe_start_coordinate <-  main_probe_start_coord 
      partner_probe_selection$main_probe_end_coordinate <- main_probe_end_coord 
      output_list[[i]] <- partner_probe_selection
    }
    
  } 
  
  output_df <- do.call("rbind", output_list) %>% dplyr::rename(partner_probe_type = probe_type, 
                                                               partner_probe_id = qseqid,
                                                               partner_probe_start_coordinate = start_coordinate,
                                                               partner_probe_end_coordinate = end_coordinate)
  
  return(output_df)
  
}


#Wrapper function to label main probes and theoretical partner probes
#input_df can be blast_res_species_99_on or blast_res_species_90_to_98_on etc
#incorporate btop results from BLAST
label_main_and_partners <- function(input_df){
  
  # Fail gracefully if input is null or empty
  if(is.null(input_df) || nrow(input_df) == 0){
    return(NULL)
  }
  
  #the union of main probes, left_partners and right partners, unique
  theoretical_partners <- generate_partner_probe_ids(input_df %>% dplyr::select(qseqid)) %>% unique(.)
  
  #Remember that these partner probes may or may not exist in probeMiner blockparse results
  #For each main probe (chosen as "main" for specificity to a target), look for existence of a valid partner probe. 
  #Note: A partner probe can also be a main probe. Main probe just refers to the probes that were in the initial shortlist
  theoretical_partners$probe_type <- ifelse(theoretical_partners$qseqid %in% input_df$qseqid,
                                            "main_probe", "partner_probe")
  
  candidate_pairs <- generate_probes_w_partners(theoretical_partners)
  
  #returns an output df with candidate probe pairs (before off-target filtering and thermoanalysis)
  return(candidate_pairs)
  
}

#function to return and summarize blast results for partner probes (adjacent to any "main probe", which were themselves chosen for higher specificity)
#input_df is the output of label_main_and_partners()
#N.B: blast_res_filt excludes all hits that have a positive scoring match of <= 20 because this means more than 5 mismatches between the 25 nt query to subject
summarize_blast_for_candidate_pairs <- function(input_df, blast_results){
  
  # Fail gracefully if input is null or empty
  if(is.null(input_df) || nrow(input_df) == 0){
    return(NULL)
  }
  
  partner_probe_ids <- input_df %>% pull(partner_probe_id)
  
  blast_res_filt_subset <- blast_results %>% dplyr::filter(qseqid %in% partner_probe_ids)
  
  blast_res_filt_list <- split(blast_res_filt_subset, f=blast_res_filt_subset$qseqid) #a named list, each name being a partner probe id
  
  blast_res_partners_summary_df <- lapply(blast_res_filt_list[partner_probe_ids], function(df){
    
    hit_species <- unique(df$subj_species)
    qseqid_uniq <- unique(df$qseqid)
    qseqid_vec <- rep(qseqid_uniq, times=length(hit_species))
    num_species_hits <- rep(length(hit_species),times=length(hit_species))
    
    #query molecule
    
    q_molecule <- sub("_[rR]ibosomal_RNA.*","", qseqid_vec)
    
    summary_out <- data.frame(q_molecule_partner=q_molecule,
                              qseqid_partner=qseqid_vec,
                              subj_species_partner=hit_species,
                              num_species_hits_partner=num_species_hits)
    return(summary_out)
    
  }) %>% do.call("rbind",.)
  
  return(blast_res_partners_summary_df)
}


#Function to extract and format query species and query genera

extract_and_format_query_taxon <- function(df){
  
  # Fail gracefully if input is null or empty
  if(is.null(df) || nrow(df) == 0){
    return(NULL)
  }
  
  output <- df
  
  output$query_species <- stringr::word(output$qseqid,1,2,sep="_")
  
  output$query_genera <- stringr::word(output$query_species,1,sep="_")
  
  output$query_species <- gsub(pattern="_", replacement=" ", output$query_species)
  
  return(output)
  
}

#Label whether the main probe is LHS or RHS

LHS_RHS_label_fn <- function(input_df) {
  
  # Fail gracefully if input is null or empty
  if(is.null(input_df) || nrow(input_df) == 0){
    return(NULL)
  }
  
  input_df %>%
    mutate(
      coord_cmp = main_probe_end_coordinate - partner_probe_end_coordinate,
      
      LHS = if_else(coord_cmp > 0, "main_probe", "partner_probe"),
      RHS = if_else(coord_cmp < 0, "main_probe", "partner_probe"),
      
      LHS_probe_id = if_else(coord_cmp > 0, main_probe_id, partner_probe_id),
      RHS_probe_id = if_else(coord_cmp < 0, main_probe_id, partner_probe_id)
    ) %>%
    dplyr::select(-coord_cmp)
  
}

#
parse_blast_qc_flags <- function(btop, aln_length, qstart, qend, qlen,
                                 side = c("LHS", "RHS"),
                                 n = 5) {
  
  side <- match.arg(side)
  
  # Coverage checks (termini/ligation junction) 
  terminal_mismatch <- FALSE
  
  # at least n=5 nucleotides (ligation junction) failed to align for RHS probe  
  if (side == "RHS" && qstart > n) {
    terminal_mismatch <- TRUE
  }
  # at least n=5 nucleotides (ligation junction) failed to align for LHS probe  
  if (side == "LHS" && (qlen - qend) >= n) {
    terminal_mismatch <- TRUE
  }
  
  # split the BTOP string into atomic tokens
  #e.g. regmatches("10AC5-T3", gregexpr("[0-9]+|[A-Z\\-][A-Z\\-]", "10AC5-T3"))[[1]] returns c("10", "AC", "5", "-T", "3")
  tokens <- regmatches(
    btop,
    gregexpr("[0-9]+|[A-Z\\-][A-Z\\-]", btop)
  )[[1]]
  
  #initialize vectors
  aln <- character()
  five_prime_unalign_vec <- rep("unalign", qstart-1) #for nts that fail to align at 5' end of probe (query)
  three_prime_unalign_vec <- rep("unalign", qlen-qend ) #for nts that fail to align at 3' end of probe (query)
  
  mismatch_count <- 0L
  
  for (tok in tokens) {
    
    # Match runs -> integers in btop
    if (grepl("^[0-9]+$", tok)) {
      aln <- c(aln, rep("match", as.integer(tok)))
      next
    }
    
    #the first character in the token is the query, the second is the subject
    q <- substr(tok, 1, 1)
    s <- substr(tok, 2, 2)
    
    # Gap in query, ignore
    if (q == "-") next
    
    # Gap in subject
    if (s == "-") {
      aln <- c(aln, "gap")
      next
    }
    
    # True mismatch because this is neither a match, a gap in query or a gap in subject
    # Appends the string "mismatch" to the aln vector
    aln <- c(aln, "mismatch")
    mismatch_count <- mismatch_count + 1L
  }
  ##############################################################################
  #combined vector with unaligned nucleotides taken into account at both ends
  combined_aln <- c(five_prime_unalign_vec, aln, three_prime_unalign_vec)
  
  unaligned_count <- qlen - aln_length
  
  mismatch_and_unaligned_count <- mismatch_count + unaligned_count
  ##############################################################################
  #Terminal mismatch/gap check if there is alignment at the ligation junction #
  #alignment too short, with fewer than n query bases aligned at all
  #Any non-match? Returns FALSE if all n bases in the terminal window are perfect matches
  if (length(aln) < n) {
    terminal_mismatch <- TRUE
  } else if (!terminal_mismatch) {
    terminal <- if (side == "RHS") head(combined_aln, n) else tail(combined_aln, n)
    if (any(terminal != "match")) {
      terminal_mismatch <- TRUE
    }
  }
  
  # Return both flags
  list(
    terminal_mismatch_or_unaligned = terminal_mismatch,
    five_or_more_mismatches_and_unaligned = mismatch_and_unaligned_count >= 5
  )
}

#
## wrapper function
# #Matches to off-target genes should have at least five mismatches in at least one of the LHS or RHS probes to prevent efficient hybridization
#probe pairs can be blast_res_species_99_on_theoretical_pairs for example
parse_LHS_and_RHS_matches <- function(full_blast_res, probe_pairs){
  
  # Fail gracefully if input is null or empty
  if(is.null(probe_pairs) || nrow(probe_pairs) == 0){
    return(NULL)
  }
  
  ########
  ##LHS###
  ########
  
  LHS_blast_res <- full_blast_res %>% dplyr::filter(qseqid %in% probe_pairs$LHS_probe_id)
  
  #Because SIMPLIFY=FALSE, LHS_qc is a list with length = nrow(LHS_blast_res), and each element is another named list
  LHS_qc <- mapply(
    parse_blast_qc_flags,
    btop   = LHS_blast_res$btop,
    aln_length = LHS_blast_res$length,
    qstart = LHS_blast_res$qstart,
    qend   = LHS_blast_res$qend,
    qlen   = LHS_blast_res$qlen,
    MoreArgs = list(side = "LHS", n = 5),
    SIMPLIFY = FALSE
  )
  
  #[[ is the list-extraction operator used as a function, "terminal_mismatch" is passed as the index e.g list[["terminal_mismatch]]
  #FUN.VALUE= logical(1) enforces each extracted value must be one logical: TRUE or FALSE
  LHS_blast_res$terminal_mismatch_or_unaligned  <- vapply(LHS_qc, `[[`, logical(1), "terminal_mismatch_or_unaligned")
  LHS_blast_res$five_or_more_mismatches_and_unaligned <- vapply(LHS_qc, `[[`, logical(1), "five_or_more_mismatches_and_unaligned")
  
  
  ########
  ##RHS###
  ########
  RHS_blast_res <- full_blast_res %>% dplyr::filter(qseqid %in% probe_pairs$RHS_probe_id)
  
  RHS_qc <- mapply(
    parse_blast_qc_flags,
    btop   = RHS_blast_res$btop,
    aln_length = RHS_blast_res$length,
    qstart = RHS_blast_res$qstart,
    qend   = RHS_blast_res$qend,
    qlen   = RHS_blast_res$qlen,
    MoreArgs = list(side = "RHS", n = 5),
    SIMPLIFY = FALSE
  )
  
  RHS_blast_res$terminal_mismatch_or_unaligned  <- vapply(RHS_qc, `[[`, logical(1), "terminal_mismatch_or_unaligned")
  RHS_blast_res$five_or_more_mismatches_and_unaligned <- vapply(RHS_qc, `[[`, logical(1), "five_or_more_mismatches_and_unaligned")
  
  
  output <- tibble::lst(LHS_blast_res, RHS_blast_res)
  
  return(output)
  
}

############################################
#helper functions for find_shared_targets()#
############################################

simplify_tax_string <- function(x, mode = c("species", "genus")) {
  # Ensures 'mode' is one of the two options
  mode <- match.arg(mode)
  
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    
    # --- Step 1: Find Genus ---
    # Look for capitalized word after the last semicolon
    tail_part <- sub(".*;", "", s)
    m <- regexpr("\\b[A-Z][a-z]+\\b", tail_part)
    
    if (m != -1) {
      genus <- regmatches(tail_part, m)
    } else {
      # Fallback: Find the LAST capitalized word in the whole string
      all_caps_list <- gregexpr("\\b[A-Z][a-z]+\\b", s)
      # FIX: regmatches returns a list; we extract the first (and only) element
      matches <- regmatches(s, all_caps_list)[[1]]
      
      if (length(matches) == 0) return(NA_character_)
      genus <- matches[length(matches)]
    }
    
    # --- Step 2: Construct Label ---
    if (mode == "genus") {
      return(genus)
    }
    
    # If mode is species:
    # Use the tail context if genus is present there, otherwise full string
    if (grepl(paste0("\\b", genus, "\\b"), tail_part)) {
      return(sub(paste0(".*\\b(", genus, "\\b.*)$"), "\\1", tail_part))
    } else {
      return(sub(paste0(".*\\b(", genus, "\\b.*)$"), "\\1", s))
    }
    
  }, character(1))
}

#count_num_terminal_no_match(df=test_90to98_on_by_probe, parsed_blast_res = blast_res_species_90to98_on_for_candidates)
count_num_terminal_no_match <- function(df, parsed_blast_res){
  
  combined_blast_res_terminal_no_match <- dplyr::bind_rows(parsed_blast_res) %>%
    dplyr::filter(terminal_mismatch_or_unaligned)
  
  #create placeholder if so
  if (nrow (combined_blast_res_terminal_no_match) == 0) {
    
    combined_blast_res_terminal_no_match <- data.frame(qseqid = "empty", terminal_mismatch_or_unaligned = TRUE)
    
  }
  
  df_off_target_species <- df %>% dplyr::filter(on_target_blast_species == FALSE)
  
  df_off_target_genus <- df %>% dplyr::filter(on_target_blast_genera == FALSE)
  
  if (nrow(df_off_target_species) > 0) {
    
    #main probe, species
    df_off_target_species_terminal_MAIN <- combined_blast_res_terminal_no_match %>% 
      dplyr::filter(qseqid %in% df_off_target_species$main_probe_id) %>%
      dplyr::select(qseqid, subj_species, terminal_mismatch_or_unaligned) %>% 
      group_by(qseqid, subj_species) %>% 
      dplyr::summarise(main_probe_terminal_mismatch_unaligned_count_PER_species=sum(terminal_mismatch_or_unaligned)) %>%
      dplyr::mutate(main_probe_terminal_mismatch_unaligned_count_unique_species = as.integer(main_probe_terminal_mismatch_unaligned_count_PER_species > 0)) %>%
      ungroup(.)
    
    
    #partner probe, species
    df_off_target_species_terminal_PARTNER <- combined_blast_res_terminal_no_match %>% 
      dplyr::filter(qseqid %in% df_off_target_species$partner_probe_id) %>%
      dplyr::select(qseqid, subj_species, terminal_mismatch_or_unaligned) %>% 
      group_by(qseqid, subj_species) %>% 
      dplyr::summarise(partner_probe_terminal_mismatch_unaligned_count_PER_species=sum(terminal_mismatch_or_unaligned)) %>%
      dplyr::mutate(partner_probe_terminal_mismatch_unaligned_count_unique_species = as.integer(partner_probe_terminal_mismatch_unaligned_count_PER_species > 0)) %>%
      ungroup(.)
    
    
    ###
    
    df_off_target_species_terminal <- left_join(x=df_off_target_species,
                                                y=df_off_target_species_terminal_MAIN,
                                                by = c("main_probe_id"="qseqid",
                                                       "shared_hits"="subj_species"))
    
    df_off_target_species_terminal <- left_join(x=df_off_target_species_terminal,
                                                y=df_off_target_species_terminal_PARTNER,
                                                by = c("partner_probe_id"="qseqid",
                                                       "shared_hits"="subj_species")) 
    
    df_off_target_species_terminal <- df_off_target_species_terminal %>%
      mutate(main_probe_terminal_mismatch_unaligned_count_PER_species = tidyr::replace_na(main_probe_terminal_mismatch_unaligned_count_PER_species, 0),
             partner_probe_terminal_mismatch_unaligned_count_PER_species = tidyr::replace_na(partner_probe_terminal_mismatch_unaligned_count_PER_species, 0),
             main_probe_terminal_mismatch_unaligned_count_unique_species = tidyr::replace_na(main_probe_terminal_mismatch_unaligned_count_unique_species,0),
             partner_probe_terminal_mismatch_unaligned_count_unique_species = tidyr::replace_na(partner_probe_terminal_mismatch_unaligned_count_unique_species,0)
      )
    
  } else if (nrow(df_off_target_species) == 0){
    
    df_off_target_species_terminal <- df_off_target_species
    df_off_target_species_terminal$main_probe_terminal_mismatch_unaligned_count_PER_species <- integer(0)    
    df_off_target_species_terminal$main_probe_terminal_mismatch_unaligned_count_unique_species <- integer(0)
    df_off_target_species_terminal$partner_probe_terminal_mismatch_unaligned_count_PER_species <- integer(0)
    df_off_target_species_terminal$partner_probe_terminal_mismatch_unaligned_count_unique_species <- integer(0)
    
  }
  
  if (nrow(df_off_target_genus) > 0) {
    
    #main probe, genus
    df_off_target_genus_terminal_MAIN <- combined_blast_res_terminal_no_match %>% 
      dplyr::filter(qseqid %in% df_off_target_genus$main_probe_id) %>%
      dplyr::select(qseqid, subj_species, terminal_mismatch_or_unaligned) %>% 
      group_by(qseqid, subj_species) %>% 
      dplyr::summarise(main_probe_terminal_mismatch_unaligned_count_PER_genus=sum(terminal_mismatch_or_unaligned)) %>%
      ungroup(.)
    #partner probe, genus
    df_off_target_genus_terminal_PARTNER <- combined_blast_res_terminal_no_match %>% 
      dplyr::filter(qseqid %in% df_off_target_genus$partner_probe_id) %>%
      dplyr::select(qseqid, subj_species, terminal_mismatch_or_unaligned) %>% 
      group_by(qseqid, subj_species) %>% 
      dplyr::summarise(partner_probe_terminal_mismatch_unaligned_count_PER_genus=sum(terminal_mismatch_or_unaligned)) %>%
      ungroup(.)
    
    ###
    
    df_off_target_genus_terminal <- left_join(x=df_off_target_genus,
                                              y=df_off_target_genus_terminal_MAIN,
                                              by = c("main_probe_id"="qseqid",
                                                     "shared_hits"="subj_species"))
    
    df_off_target_genus_terminal <- left_join(x=df_off_target_genus_terminal,
                                              y=df_off_target_genus_terminal_PARTNER,
                                              by = c("partner_probe_id"="qseqid",
                                                     "shared_hits"="subj_species")) 
    
    df_off_target_genus_terminal <- df_off_target_genus_terminal %>%
      mutate(main_probe_terminal_mismatch_unaligned_count_PER_genus = tidyr::replace_na(main_probe_terminal_mismatch_unaligned_count_PER_genus, 0),
             partner_probe_terminal_mismatch_unaligned_count_PER_genus = tidyr::replace_na(partner_probe_terminal_mismatch_unaligned_count_PER_genus, 0))
    
  } else if (nrow(df_off_target_genus) == 0){
    
    df_off_target_genus_terminal <- df_off_target_genus
    
  }
  
  #Generate summaries
  
  
  if (nrow(df_off_target_species_terminal) > 0) {
    
    df_off_target_species_terminal_summary <- df_off_target_species_terminal %>% 
      dplyr::select(main_probe_id, partner_probe_id,
                    main_probe_terminal_mismatch_unaligned_count_PER_species, 
                    partner_probe_terminal_mismatch_unaligned_count_PER_species,
                    main_probe_terminal_mismatch_unaligned_count_unique_species,
                    partner_probe_terminal_mismatch_unaligned_count_unique_species) %>%
      group_by(main_probe_id, partner_probe_id) %>% 
      summarise(main_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target=sum(main_probe_terminal_mismatch_unaligned_count_PER_species),
                partner_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target=sum(partner_probe_terminal_mismatch_unaligned_count_PER_species),
                main_probe_terminal_mismatch_unaligned_num_uniq_species_off_target=sum(main_probe_terminal_mismatch_unaligned_count_unique_species),
                partner_probe_terminal_mismatch_unaligned_num_uniq_species_off_target=sum(partner_probe_terminal_mismatch_unaligned_count_unique_species)
      ) %>% 
      ungroup(.)
    
  } else if (nrow(df_off_target_species_terminal) == 0){
    
    df_off_target_species_terminal_summary <- data.frame(main_probe_id = "empty",
                                                         partner_probe_id = "empty",
                                                         main_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target=0,
                                                         partner_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target=0,
                                                         main_probe_terminal_mismatch_unaligned_num_uniq_species_off_target=0,
                                                         partner_probe_terminal_mismatch_unaligned_num_uniq_species_off_target=0)
    
  }
  
  if (nrow(df_off_target_genus_terminal) > 0) {
    
    df_off_target_genus_terminal_summary <- df_off_target_genus_terminal %>% 
      dplyr::select(main_probe_id, partner_probe_id,
                    main_probe_terminal_mismatch_unaligned_count_PER_genus, 
                    partner_probe_terminal_mismatch_unaligned_count_PER_genus) %>%
      group_by(main_probe_id, partner_probe_id) %>% 
      summarise(main_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target=sum(main_probe_terminal_mismatch_unaligned_count_PER_genus),
                partner_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target=sum(partner_probe_terminal_mismatch_unaligned_count_PER_genus)) %>% 
      ungroup(.)
    
  } else if (nrow(df_off_target_genus_terminal) == 0){
    
    df_off_target_genus_terminal_summary <- data.frame(main_probe_id = "empty",
                                                       partner_probe_id = "empty",
                                                       main_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target=0,
                                                       partner_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target=0)
    
  }
  
  
  output_summary <- full_join(df_off_target_species_terminal_summary, df_off_target_genus_terminal_summary,
                              by = c("main_probe_id", "partner_probe_id")) %>% 
    dplyr::filter(main_probe_id != "empty") %>%  replace(is.na(.), 0)
  
  return(output_summary) #will return an empty dataframe if there are no terminal mismatches or terminal unaligned nucleotides in the blast results
  
}


################################################
########Find shared hits for probe pairs  ######
################################################

#for each main probe with its partner, find the common on-targets and off-targets (shared with the partner) and count the number of species that will be hit, for both on targets and off targets

#this function is downstream of parse_LHS_and_RHS_matches(), which outputs a list object for blast_res_list
#probe_pairs can be blast_res_species_99_on_theoretical_pairs etc.
#candidate_blast_results can be blast_res_species_99_on_for_candidates or equivalent

#We want to find both shared targets at the species level and also at the genome level, because some genomes may be mis-annotated

find_shared_targets <- function(blast_res_list, probe_pairs) {
  
  # Fail gracefully if input is null or empty
  if(is.null(blast_res_list) || is.null(probe_pairs) || nrow(probe_pairs) == 0){
    return(list(paired_hits_by_probe = NULL, paired_genome_hits_by_probe = NULL, probe_pair_hit_summary = NULL, probe_pair_hit_summary_fmt = NULL))
  }
  
  #Analyze LHS BLAST results, filtering for "hits", which can be on target or off target
  LHS_blast_res <- blast_res_list[["LHS_blast_res"]] %>% dplyr::filter(five_or_more_mismatches_and_unaligned == FALSE)
  #Analyze RHS BLAST results
  RHS_blast_res <- blast_res_list[["RHS_blast_res"]] %>% dplyr::filter(five_or_more_mismatches_and_unaligned == FALSE)
  
  #########################
  #Shared species analysis#
  #########################
  # Pre-split BLAST hits by probe for fast lookup
  species_by_LHS_probe <- LHS_blast_res %>%
    distinct(qseqid, subj_species) %>%
    group_by(qseqid) %>%
    summarise(species = list(subj_species), .groups = "drop")
  
  species_by_RHS_probe <- RHS_blast_res %>%
    distinct(qseqid, subj_species) %>%
    group_by(qseqid) %>%
    summarise(species = list(subj_species), .groups = "drop")
  
  
  # Join species lists to probe pairs
  paired_hits_by_probe <- probe_pairs %>% 
    inner_join(species_by_LHS_probe, by = c("LHS_probe_id" = "qseqid")) %>%
    rename(LHS_hits = species) %>%
    inner_join(species_by_RHS_probe, by = c("RHS_probe_id" = "qseqid")) %>%
    rename(RHS_hits = species) %>%
    mutate(shared_hits = Map(intersect, LHS_hits, RHS_hits)) %>%
    select(main_probe_id, partner_probe_id, shared_hits, LHS, RHS) %>%
    unnest(shared_hits) %>% unique(.)
  
  
  paired_hits_by_probe$subj_species <-  simplify_tax_string(paired_hits_by_probe$shared_hits, mode = "species")
  
  paired_hits_by_probe$subj_species <- gsub(pattern=" ", replacement="_", x = paired_hits_by_probe$subj_species)
  
  paired_hits_by_probe$subj_genera <-  simplify_tax_string(paired_hits_by_probe$shared_hits, mode = "genus")
  
  ###############################################
  #detect on and off targets for shared species##
  ###############################################
  
  #to avoid modifying the actual main_probe_id, create a temp column
  paired_hits_by_probe$main_probe_id_temp <- str_replace_all(paired_hits_by_probe$main_probe_id, 
                                                             species_rename_map_underscore)
  
  paired_hits_by_probe$subj_genera <- if_else(str_detect(paired_hits_by_probe$subj_species, Phocaeicola_mask),
                                              "Phocaeicola",
                                              paired_hits_by_probe$subj_genera)
  
  paired_hits_by_probe$subj_genera <- if_else(str_detect(paired_hits_by_probe$subj_species, Cutibacterium_mask),
                                              "Cutibacterium",
                                              paired_hits_by_probe$subj_genera)
  
  
  paired_hits_by_probe <- paired_hits_by_probe %>%
    mutate(subj_species = if_else(
      str_detect(subj_species, Phocaeicola_mask),
      str_replace(subj_species, pattern="Bacteroides", replacement="Phocaeicola"),
      subj_species
    )) %>% mutate(subj_species = if_else(
      str_detect(subj_species, Cutibacterium_mask),
      str_replace(subj_species, pattern="Propionibacterium", replacement="Cutibacterium"),
      subj_species
    )) 
  
  
  #remove underscores to facilitate comparison
  paired_hits_by_probe$subj_species <- gsub(pattern="_", replacement=" ", x = paired_hits_by_probe$subj_species)

  #Use stringr to extract the species/genera
  paired_hits_by_probe$query_species <- stringr::word(paired_hits_by_probe$main_probe_id,1,2,sep="_")

  #Replace underscores
  paired_hits_by_probe$query_species <- gsub(pattern="_", replacement=" ",paired_hits_by_probe$query_species)
  
  paired_hits_by_probe$on_target_blast_species <- ifelse(str_detect(string=paired_hits_by_probe$subj_species, 
                                                                    pattern=paired_hits_by_probe$query_species), TRUE, FALSE)
  
  paired_hits_by_probe$on_target_blast_genera <- ifelse(str_detect(string=paired_hits_by_probe$main_probe_id_temp, 
                                                                   pattern=paired_hits_by_probe$subj_genera), TRUE, FALSE)
  
  #drop the temp column
  paired_hits_by_probe <- paired_hits_by_probe %>% dplyr::select(-main_probe_id_temp)
  
  
  #########################
  #Shared genomes analysis#
  #########################
  # Pre-split BLAST hits by probe for fast lookup
  genomes_by_LHS_probe <- LHS_blast_res %>%
    distinct(qseqid, stitle) %>%
    group_by(qseqid) %>%
    summarise(subj_genome = list(stitle), .groups = "drop")
  
  genomes_by_RHS_probe <- RHS_blast_res %>%
    distinct(qseqid, stitle) %>%
    group_by(qseqid) %>%
    summarise(subj_genome = list(stitle), .groups = "drop")
  
  # Join genome lists to probe pairs
  paired_genome_hits_by_probe <- probe_pairs %>% 
    inner_join(genomes_by_LHS_probe, by = c("LHS_probe_id" = "qseqid")) %>%
    dplyr::rename(LHS_genome_hits = subj_genome) %>%
    inner_join(genomes_by_RHS_probe, by = c("RHS_probe_id" = "qseqid")) %>%
    dplyr::rename(RHS_genome_hits = subj_genome) %>%
    mutate(shared_genome_hits = Map(intersect, LHS_genome_hits, RHS_genome_hits)) %>%
    select(main_probe_id, partner_probe_id, shared_genome_hits, LHS, RHS) %>%
    unnest(shared_genome_hits) %>% unique(.)
  
  
  ####################################################
  ##SUMMARY STATS###
  
  #Return the formatted BLAST results and a hit summary of 
  #num_on_target_species, num_off_target_species, num_on_target_genera, num_off_target_genera for each probe pair
  
  probe_pair_hit_summary <- paired_hits_by_probe  %>% 
    group_by(main_probe_id, partner_probe_id, LHS, RHS) %>%
    summarise(num_on_target_species=sum(on_target_blast_species==TRUE),
              num_off_target_species=sum(on_target_blast_species==FALSE),
              num_on_target_genera=sum(on_target_blast_genera==TRUE),
              num_off_target_genera=sum(on_target_blast_genera==FALSE)) %>% ungroup(.)
  
  #add number of off target genomes (IMPORTANT DEFINITION) for which there are mismatched or unaligned bases in the termini (ligation junction), for both species and genus level matches
  #for each off target, count the number of instances where there is a mismatch or unaligned base at the ligation junction (+- 5 nt) separately for the main and partner probe
  
  probe_pair_off_target_termini_summary <- count_num_terminal_no_match(df=paired_hits_by_probe, 
                                                                       parsed_blast_res = blast_res_list)
  
  if (nrow(probe_pair_off_target_termini_summary > 0)){
    
    probe_pair_hit_summary_fmt <- left_join(probe_pair_hit_summary, 
                                            probe_pair_off_target_termini_summary,
                                            by = c("main_probe_id", "partner_probe_id")) %>%  replace(is.na(.), 0)
    
  } else {
    
    probe_pair_hit_summary_fmt <- probe_pair_hit_summary
    probe_pair_hit_summary_fmt$main_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target <- 0
    probe_pair_hit_summary_fmt$partner_probe_terminal_mismatch_unaligned_count_species_GENOMES_off_target <- 0
    probe_pair_hit_summary_fmt$main_probe_terminal_mismatch_unaligned_num_uniq_species_off_target <- 0
    probe_pair_hit_summary_fmt$partner_probe_terminal_mismatch_unaligned_num_uniq_species_off_target <- 0
    probe_pair_hit_summary_fmt$main_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target <- 0
    probe_pair_hit_summary_fmt$partner_probe_terminal_mismatch_unaligned_count_genus_GENOMES_off_target <- 0
  }
  
  
  output <- tibble::lst(paired_hits_by_probe, paired_genome_hits_by_probe, probe_pair_hit_summary_fmt)
  
  return(output)
  
}

#filter probes (blast_res_species_99_on_candidate_pairs$probe_pair_hit_summary_fmt) for only those that have passed thermoanalysis

filter_by_thermoanalysis <- function(df){
  
  # Handle NULL or empty input
  if (is.null(df) || nrow(df) == 0 || !all(c("main_probe_id", "partner_probe_id") %in% names(df))) {
    return(NULL)
  }
  
  df$main_probe_thermo_pass <- ifelse(df$main_probe_id %in% probe_primer3_passed$Sequence_ID, TRUE, FALSE)
  df$partner_probe_thermo_pass <- ifelse(df$partner_probe_id %in% probe_primer3_passed$Sequence_ID, TRUE, FALSE)
  
  df_filt <- df %>% dplyr::filter(main_probe_thermo_pass==TRUE & partner_probe_thermo_pass==TRUE)
  
  # Return NULL if no probes pass filters
  if (nrow(df_filt) == 0) {
    return(NULL)
  }
  
  #add the actual probe target sequences into the summary file
  #reminder: probe target sequences =/= probe sequences
  #probably should add the Tm values here...
  probe_sequences <- probe_primer3_passed %>% dplyr::select(Sequence_ID, Sequence)
  
  df_filt <- inner_join(df_filt, probe_sequences, by = c("main_probe_id" = "Sequence_ID")) %>% 
    dplyr::rename(main_probe_target_seq=Sequence) %>% 
    inner_join(., probe_sequences, by = c("partner_probe_id" = "Sequence_ID")) %>%
    dplyr::rename(partner_probe_target_seq=Sequence)  
  
  return(df_filt)
  
}

revcomp <- function(seqs) {
  stopifnot(is.character(seqs))
  
  # Mapping table
  comp <- c(
    "A"="T", "T"="A",
    "C"="G", "G"="C",
    "a"="T", "t"="A",
    "c"="G", "g"="C"
  )
  
  vapply(seqs, function(s) {
    bases <- strsplit(s, "", fixed = TRUE)[[1]]
    rc <- rev(comp[bases])
    paste(rc, collapse = "")
  }, character(1))
}




#function to reverse complement and add Visium adapters

append_Visium_adapters <- function(df_fmt){
  
  #5' end LHS adapter
  LHS_adapter <- "CCTTGGCACCCGAGAATTCCA"
  
  #3' end RHS tail with 30 As
  RHS_tail <- rep("A", times = 30) %>% paste(., sep="", collapse="") 
  
  #the actual reverse complement
  df_fmt$main_probe_seq <- revcomp(df_fmt$main_probe_target_seq)
  df_fmt$partner_probe_seq <- revcomp(df_fmt$partner_probe_target_seq)
  
  #extract LHS and RHS probes into separate dataframes
  LHS_main_probe <- df_fmt %>% 
    dplyr::filter(LHS=="main_probe") %>% 
    dplyr::select(main_probe_seq, main_probe_id) %>%
    dplyr::rename(LHS_probe_sequence= main_probe_seq,
                  LHS_probe_id= main_probe_id)
  
  LHS_partner_probe <- df_fmt %>% 
    dplyr::filter(LHS=="partner_probe") %>% 
    dplyr::select(partner_probe_seq, partner_probe_id) %>%
    dplyr::rename(LHS_probe_sequence= partner_probe_seq,
                  LHS_probe_id= partner_probe_id)
  
  RHS_main_probe <- df_fmt %>% 
    dplyr::filter(RHS=="main_probe") %>% 
    dplyr::select(main_probe_seq, main_probe_id) %>%
    dplyr::rename(RHS_probe_sequence= main_probe_seq,
                  RHS_probe_id= main_probe_id)
  
  RHS_partner_probe <- df_fmt %>%
    dplyr::filter(RHS=="partner_probe") %>% 
    dplyr::select(partner_probe_seq, partner_probe_id) %>%
    dplyr::rename(RHS_probe_sequence= partner_probe_seq,
                  RHS_probe_id= partner_probe_id)
  
  LHS_all_probes <- rbind(LHS_main_probe, LHS_partner_probe)
  RHS_all_probes <- rbind(RHS_main_probe, RHS_partner_probe)
  
  #add visium adapters
  LHS_all_probes$LHS_probe_sequence <- paste0(LHS_adapter,LHS_all_probes$LHS_probe_sequence)
  
  RHS_all_probes$RHS_probe_sequence <- paste0(RHS_all_probes$RHS_probe_sequence,RHS_tail)
  
  output <- tibble::lst(LHS_all_probes, RHS_all_probes)
  
  return(output)
  
}


#Probes overlap if y1 <= x2, assuming the probed regions are arranged in order AND of equal size
#https://stackoverflow.com/questions/41747742/collapse-rows-with-overlapping-ranges

#this function is per species
select_species_probes_from_shortlist <- function(df, max_probes = 7) {
  
  # 1) Compute genomic intervals
  df <- df %>%
    mutate(
      range_min = pmin(main_probe_start_coordinate,
                       main_probe_end_coordinate,
                       partner_probe_start_coordinate,
                       partner_probe_end_coordinate),
      range_max = pmax(main_probe_start_coordinate,
                       main_probe_end_coordinate,
                       partner_probe_start_coordinate,
                       partner_probe_end_coordinate)
    )
  
  out <- df %>%
    group_by(query_molecule) %>%
    group_modify(~ {
      x <- .x %>% arrange(range_min)
      
      if (nrow(x) == 0) return(x)
      
      # 2) Assign overlap groups
      x <- x %>%
        mutate(
          probe_group = cumsum(
            lag(range_max, default = first(range_max)) < range_min
          )
        )
      
      # 3) Pick best probe per overlapping region
      best_per_group <- x %>%
        group_by(probe_group) %>%
        arrange(num_off_target_species) %>%
        slice(1) %>%
        ungroup()
      
      # 4) Remove overlaps between selected probes
      best_nonoverlap <- best_per_group %>%
        arrange(range_min) %>%
        filter(range_min >= cummax(lag(range_max, default = -Inf)))
      
      # 5) Pick top N best probes
      selected <- best_nonoverlap %>%
        arrange(num_off_target_species) %>%
        slice_head(n = max_probes) %>%
        mutate(shortlisted = TRUE,
               non_overlapping_region = TRUE)
      
      # 6) Mark rejected probes
      rejected <- x %>%
        anti_join(selected, by = "partner_probe_id") %>%
        mutate(shortlisted = FALSE,
               non_overlapping_region = FALSE)
      
      bind_rows(selected, rejected)
    }) %>%
    ungroup()
  
  out
}



#this function for shortlisting genus probes
select_genus_probes_from_shortlist <- function(df, max_probes = 7) {
  
  # 1) Compute interval
  df <- df %>%
    mutate(
      range_min = pmin(main_probe_start_coordinate,
                       main_probe_end_coordinate,
                       partner_probe_start_coordinate,
                       partner_probe_end_coordinate),
      range_max = pmax(main_probe_start_coordinate,
                       main_probe_end_coordinate,
                       partner_probe_start_coordinate,
                       partner_probe_end_coordinate)
    )
  
  # 2) Process each molecule independently
  out <- df %>%
    group_by(query_molecule) %>%
    group_modify(~ {
      x <- .x %>% arrange(range_min)
      
      if (nrow(x) == 0) return(x)
      
      # 3) Assign overlap groups
      x <- x %>%
        mutate(
          probe_group = cumsum(
            lag(range_max, default = first(range_max)) < range_min
          )
        )
      
      # 4) Pick best probe per overlapping group
      best_per_group <- x %>%
        group_by(probe_group) %>%
        arrange(num_off_target_genera,
                desc(num_on_target_genera)) %>%
        slice(1) %>%
        ungroup()
      
      # 5) Remove overlaps between chosen probes
      best_nonoverlap <- best_per_group %>%
        arrange(range_min) %>%
        filter(range_min >= cummax(lag(range_max, default = -Inf)))
      
      # 6) Pick top N
      selected <- best_nonoverlap %>%
        arrange(num_off_target_genera,
                desc(num_on_target_genera)) %>%
        slice_head(n = max_probes) %>%
        mutate(shortlisted = TRUE,
               non_overlapping_region = TRUE)
      
      # 7) Mark rejected
      rejected <- x %>%
        anti_join(selected, by = "partner_probe_id") %>%
        mutate(shortlisted = FALSE,
               non_overlapping_region = FALSE)
      
      bind_rows(selected, rejected)
    }) %>%
    ungroup()
  
  out
}
