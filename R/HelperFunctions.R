#################################################
###############Extract Feature Tags##############
#################################################
#' Internal Function To Extract Annotation Features
#' @param text Character. Chunk of GenBank annotations used for searching
#' @param tag Character. Tag to search for.
#' @param remove_spaces Logical. Should spaces be removed?
#' @noRd
#' @author Samuel R. Borstein

extract_gb_tag <- function(text, tag, remove_spaces = FALSE) {
  pattern <- paste0("/", tag, "=(\"[^\"]+\"|[^\\s]+)")#general search pattern
  extracted <- stringr::str_extract(text, pattern)#Find relevant line and isolate
  if (is.na(extracted)) return(NA)#If not extracted, make NA
  cleaned <- stringr::str_replace(extracted, paste0("/", tag, "="), "")#remove annotation tag
  cleaned <- stringr::str_replace_all(cleaned, "\"", "")#remove non-target text
  if (remove_spaces == TRUE){#If to remove spaces control flow. We choose not to by default and hard code it.
    cleaned <- gsub("\\s+", "", cleaned)#if there are spaces, collapse them.
  }
  return(cleaned)
}

############################################################
################Parse Coordinates & Joins###################
############################################################
#' Internal Parser Function
#' @param coord_str Character vector representing a coordinate string of a GenBank accession. Obtained via parse_genbank.
#' @param AccNo description
#' @examples
#' Accessions <- c("PV747861.1")
#' raw_gb_text <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "gbwithparts", retmode = "text")#Use entrez to access accession
#' gb_lines <- strsplit(raw_gb_text, "\n|\r\n")[[1]]#string split lines
#' parse_genbank(ReadGB = gb_lines)
#' @noRd
#' @author Samuel R. Borstein
parse_gb_coordinates <- function(coord_str, primary_accession) {
  coord_str <- gsub("\\s+", "", coord_str)#Remove spacing between ranges
  is_global_comp <- grepl("^complement\\(join\\(", coord_str)#Check for complements
  
  inner_str <- coord_str#Set coord_str as innr_str
  if (is_global_comp) {#if a global complement
    inner_str <- sub("^complement\\(join\\(", "", inner_str)#Remove leading text
    inner_str <- sub("\\)\\)$", "", inner_str)#remove anchor )
  } else if (grepl("^(join|order)\\(", inner_str)) {#if not a global complement
    inner_str <- sub("^(join|order)\\(", "", inner_str)#remove text around
    inner_str <- sub("\\)$", "", inner_str)#remove anchor )
  }
  
  exon_strs <- strsplit(inner_str, ",")[[1]]#String split to chunks
  if (is_global_comp) exon_strs <- rev(exon_strs)#if global complement
  
  parsed_exons <- list()#Initialize storage for parsed exons
  for (j in seq_along(exon_strs)) {#for each potential exon
    e_str <- exon_strs[j]#get positions
    e_strand <- ifelse(is_global_comp || grepl("complement", e_str), "-", "+")#ifelse for strand determination based on presence of complement or global complement
    
    e_seqname <- primary_accession#set sequence name to be Accession Number
    if (grepl(":", e_str)) {#if colon is found
      e_seqname <- sub(":.*", "", e_str)##remove ... for getting positions
      e_seqname <- gsub("[^A-Za-z0-9_.]", "", e_seqname)#make sure Accession is just accession
    }
    
    nums <- as.numeric(unlist(regmatches(e_str, gregexpr("\\d+", e_str))))#force numeric for numbers and getting positions
    if (length(nums) >= 1) {#if 1 or more generate a list with the name, start position, end position, strand +/-, intron/exon number, and if a join 
      parsed_exons[[length(parsed_exons) + 1]] <- list(
        seqnames = e_seqname,
        start = min(nums),
        end = max(nums),
        strand = e_strand,
        number = j,  
        loctype = ifelse(grepl("join", coord_str), "join", "normal")
      )
    }
  }
  return(dplyr::bind_rows(parsed_exons))#return the 
}

###############################################################
##############Main Internal Parser Function####################
###############################################################

#' Internal Parser Function
#' @param ReadGB Character vector of read in GenBank accession.
#' @examples
#' raw_gb_text <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "gbwithparts", retmode = "text")#Use entrez to access accession
#' gb_lines <- strsplit(raw_gb_text, "\n|\r\n")[[1]]#string split lines
#' parse_genbank(ReadGB = gb_lines)
#' @noRd
#' @author Samuel R. Borstein
#' @importFrom magrittr %>%

parse_genbank <- function(ReadGB, primary_accession) {
  orig_start <- grep("^ORIGIN", ReadGB)#Find origin flag for the sequence
  end_idx <- grep("^//", ReadGB)#Find end of sequence flag
  if (length(orig_start) > 0 && length(end_idx) > 0) {#if we have these found
    #remove line numbers and past into a single string from origin to end
    clean_dna <- gsub("[0-9 ]", "", paste(ReadGB[(orig_start + 1):(end_idx[length(end_idx)] - 1)], collapse = ""))
    full_sequence <- Biostrings::DNAStringSet(clean_dna)#Make this a DNAStringSet object
    names(full_sequence) <- primary_accession#Set accession to be the name
  } else stop("Could not find ORIGIN sequence.")#Kill if origin not found
  
  
  feat_start <- grep("^FEATURES", ReadGB)#Find the start of the features
  feat_lines <- ReadGB[(feat_start + 1):(orig_start - 1)]#Add one for start of Annotation location
  #Need to allow numbers, hyphens, to identify feature boundaries
  new_feat_idx <- grep("^ {5}[A-Za-z0-9_'-]+", feat_lines)##Match 5 spaces and alphabetical to defining features

  all_rows <- list()#Generate Empty lst
  
  #Extract basic features
  for (i in seq_along(new_feat_idx)) {#For each feature
    start_idx <- new_feat_idx[i]#Grab starting position
    end_idx <- ifelse(i < length(new_feat_idx), new_feat_idx[i+1] - 1, length(feat_lines))#Grab the end of current feature based on start of next feature
    chunk <- feat_lines[start_idx:end_idx]#subset the current feature out
    
    feat_type <- stringr::str_trim(stringr::str_extract(chunk[1], "^ {5}[A-Za-z0-9_'-]+"))#Grab current feature type from annotation
    chunk_collapsed <- paste(stringr::str_trim(chunk), collapse = " ")#Collapse feature into single line
    
    #isolate coordinates
    coord_str <- sub("^[A-Za-z0-9_'-]+\\s+", "", chunk_collapsed)#Sub out type prior to info
    coord_str <- stringr::str_trim(sub("/.*", "", coord_str))#Sub out spaces
    
    gene <- extract_gb_tag(chunk_collapsed, "gene")#Extract gene name
    locus_tag <- extract_gb_tag(chunk_collapsed, "locus_tag")#Extract Locus tag
    product <- extract_gb_tag(chunk_collapsed, "product")#Extract Product
    protein_id <- extract_gb_tag(chunk_collapsed, "protein_id")#Extract protein ID
    transcript_id_tag <- extract_gb_tag(chunk_collapsed, "transcript_id") 
    translation <- extract_gb_tag(chunk_collapsed, "translation", remove_spaces = TRUE)#If CDS, get the translation
    note <- extract_gb_tag(chunk_collapsed, "note")#Capture any notes
    
    coords_df <- parse_gb_coordinates(coord_str, primary_accession)#parse coordinates
    
    if (nrow(coords_df) > 0) {#populate the coordinate tibble for each feature type with inde, gene, tag, product, protein ID, transcipt, translation, and any relevant notes
      coords_df$feat_index <- i  
      coords_df$type <- feat_type
      coords_df$gene <- gene
      coords_df$locus_tag <- locus_tag
      coords_df$product <- product
      coords_df$protein_id <- protein_id
      coords_df$transcript_id_tag <- transcript_id_tag 
      coords_df$translation <- translation
      coords_df$note <- note
      
      all_rows[[length(all_rows) + 1]] <- coords_df
    }
  }
  
  master_df <- suppressWarnings(dplyr::bind_rows(all_rows))#Merge rows together for feature types
  
  #Check Columns to guarantee all feature tracking columns exist (prevents 'object not found' errors)
  expected_cols <- c("gene", "locus_tag", "protein_id", "transcript_id_tag")
  for (col in expected_cols) {
    if (!(col %in% colnames(master_df))) {
      master_df[[col]] <- NA
    }
  }
  
  #Strip explicit exons as they ruin trans-splicing logic so they can be inferred later.
  master_df <- master_df %>% dplyr::filter(type != "exon")
  
  #Assign genbankr like transcript IDs for keeping track internally
  master_df <- master_df %>%
    dplyr::mutate(
      #Prioritize explicit eukaryotic isoform IDs prior to falling back to locus_tag and gene
      base_name = dplyr::coalesce(transcript_id_tag, protein_id, locus_tag, gene, paste0("unnamed_", type, "_", feat_index))
    ) %>%
    dplyr::group_by(base_name) %>%
    dplyr::arrange(feat_index) %>%
    dplyr::mutate(
      instance = cumsum(c(1, diff(feat_index) > 15))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      transcript_id = dplyr::if_else(
        !is.na(transcript_id_tag) | !is.na(protein_id), 
        base_name,
        dplyr::if_else(
          !is.na(gene) & !is.na(locus_tag), 
          paste0(gene, ".", instance, "_", locus_tag), 
          paste0(base_name, ".", instance)
        )
      )
    )
  
  #Protect biological joins, sequentially number independent pieces
  master_df <- master_df %>%
    dplyr::group_by(transcript_id, type) %>%
    dplyr::arrange(feat_index, number, .by_group = TRUE) %>%
    dplyr::mutate(
      number = if (any(loctype == "join")) number else as.integer(dplyr::row_number()),
      part_id = paste0(transcript_id, "_", type, "_", number)
    ) %>%
    dplyr::ungroup()
  
  #Exon inference and tracking
  cds_rna_df <- master_df %>% dplyr::filter(type %in% c("CDS", "mRNA", "tRNA", "rRNA"))
  if (nrow(cds_rna_df) > 0) {
    inferred_exons <- cds_rna_df
    inferred_exons$type <- "exon"
    inferred_exons$part_id <- paste0(inferred_exons$transcript_id, "_exon_", inferred_exons$number)
    master_df <- dplyr::bind_rows(master_df, inferred_exons)
  }
  
  #mRNA inference
  if (!("mRNA" %in% master_df$type)) {
    cds_df <- master_df %>% dplyr::filter(type == "CDS")
    if (nrow(cds_df) > 0) {
      inferred_mrna <- cds_df
      inferred_mrna$type <- "mRNA"
      inferred_mrna$part_id <- sub("_CDS_", "_mRNA_", inferred_mrna$part_id)
      master_df <- dplyr::bind_rows(master_df, inferred_mrna)
    }
  }
  
  #Convert to standard GRanges
  master_df <- as.data.frame(master_df)
  master_gr <- GenomicRanges::makeGRangesFromDataFrame(master_df, keep.extra.columns = TRUE)
  main_types <- c("gene", "CDS", "exon", "mRNA")
  
  #Build the core output list
  final_output <- list(
    sequence = full_sequence,
    gene = master_gr[master_gr$type == "gene"],
    CDS = master_gr[master_gr$type == "CDS"],
    exon = master_gr[master_gr$type == "exon"],
    mRNA = master_gr[master_gr$type == "mRNA"]
  )
  
  #Dynamically extract and append any other feature types to list
  other_gr <- master_gr[!master_gr$type %in% main_types]
  if (length(other_gr) > 0) {
    other_split <- as.list(GenomicRanges::split(other_gr, other_gr$type))
    final_output <- c(final_output, other_split)
  }
  
  return(final_output)#Return final table of features
}

######################################
#####RANGE SEQ EXTRACT FUNCTION#######
######################################
#' Extract Sequences Based on Parsed Ranges
#' @param full_sequence The full sequence object as a DNAStringSet.
#' @param target_gr A GRanges object containing coordinate ranges for extraction.
#' @examples
#' raw_gb_text <- rentrez::entrez_fetch(db = "nuccore", id = accession, rettype = "gbwithparts", retmode = "text")#Use entrez to access accession
#' gb_lines <- strsplit(raw_gb_text, "\n|\r\n")[[1]]#string split lines
#' parse_genbank(ReadGB = gb_lines)
#' @noRd
#' @author Samuel R. Borstein
#' @import BSgenome
extract_ranges_seq <- function(full_sequence, target_gr) {
  
  #Check there is a transcript_id column for splitting joined pieces
  if (!("transcript_id" %in% colnames(GenomicRanges::mcols(target_gr)))) {
    #If provided as a raw/custom GRanges without our standard columns, 
    #just number each row individually so it doesn't crash.
    target_gr$transcript_id <- paste0("custom_range_", seq_along(target_gr))
  }
  
  #Split into a GRangesList by transcript_id
  #Needed to keep multi-exons together while keeping things distinct
  gr_list <- GenomicRanges::split(target_gr, target_gr$transcript_id)
  
  #Extract and/or stitch target sequence(s)
  seqs <- GenomicFeatures::extractTranscriptSeqs(x = full_sequence, transcripts = gr_list)
  
  return(seqs)#return the sequence
}