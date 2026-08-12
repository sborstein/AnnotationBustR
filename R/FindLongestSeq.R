#' Find the longest sequence for each species from a list of GenBank accession numbers.
#' @param Accessions A vector of GenBank accession numbers.
#' @param BatchSize Numeric. If the number of accessions is over the number provided, requests will be sent in batches of this amount. This is necessary for the NCBI servers. If you receive an HTTP 414 error, try to reduce the size of the batch. Default is 300.
#' @details For a set of GenBank accession numbers, this will return the longest sequence for in the set for species.
#' @return A list of genbank accessions numbers for the longest sequence for each taxon in a list of accession numbers.
#' @examples 
#' #a vector of 4 genbank accessions, there are two accessions for each species.
#' genbank.accessions<-c("KP978059.1","KP978060.1","JX516105.1","JX516111.1")
#' @examplesIf interactive()
#' #returns the longest sequence respectively for the two species.
#' long.seq.result <- FindLongestSeq(genbank.accessions)
#' @export

FindLongestSeq <- function(Accessions, BatchSize = 300){
  BatchAccessions <- length(Accessions)#Calc the length of the accessions for batching checks
  if(BatchAccessions > BatchSize){#If accessions is more than the batch size, batch
    BatchRequests <- split(Accessions, ceiling(seq_along(Accessions)/BatchSize))#split into batches by batch size
  }else{
    BatchRequests <- list(Accessions)#single batch if less than BatchSize threshold, but as list to reduce code to work on one or more batches
  }
  
  BatchResults <- list()#Make an empty list to stor the batches
  
  for (BatchIndex in 1:length(BatchRequests)) {
    #Obtain summaries for the current batch while checking for errors
    BatchSummary <- tryCatch({
      rentrez::entrez_summary(db = "nuccore", id = BatchRequests[[BatchIndex]])
    }, error = function(e) {
      message(paste("Error in batch", BatchIndex, ":", e$message))
      return(NULL)
    }
    )
    if (!is.null(BatchSummary)) {#If BatchSummary is non-null
      batch_df <- data.frame(
        Species = gsub(pattern = " ", replacement = "_", rentrez::extract_from_esummary(BatchSummary, "organism")),#extract organism name
        Accession = rentrez::extract_from_esummary(BatchSummary, "caption"),#Extract Accession Numbers
        Length = rentrez::extract_from_esummary(BatchSummary, "slen"),#Extract sequence langth
        stringsAsFactors = FALSE
      )
      BatchResults[[BatchIndex]] <- batch_df# Store the batch in list

    }
    Sys.sleep(0.5)#Sleep between batches to not violate NCBI terms and conditions 
  }
  FinalSequences <- do.call(rbind, BatchResults)#Bind batches together.
  
  uni.taxa <- unique(FinalSequences$Species)#get unique taxa
  long.seqs <- data.frame(matrix(nrow=length(uni.taxa), ncol=dim(FinalSequences)[2]))#empty data frame to store results
  colnames(long.seqs) <- colnames(FinalSequences)#make column names to match
  
  for (taxa.index in 1:length(uni.taxa)){
    current.tax <- subset(FinalSequences, FinalSequences$Species==uni.taxa[taxa.index])
    longest.seq <- subset(current.tax,current.tax$Length==sort(as.numeric(current.tax$Length),decreasing = TRUE)[1])[1,]#id and grab longest seq
    long.seqs[taxa.index, ] <- longest.seq#append data frame with longest seq
  }
  return(long.seqs)#Return final dataframe with longest sequences
}