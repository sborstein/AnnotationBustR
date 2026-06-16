#' Merging of two tables containing search terms to expand search term database for the AnnotationBust function.
#' @param ... the data frames of search terms you want to combine into a single data frame The Data frame(s) should have stringsAsFactors=FALSE listed if you want to sort them.
#' @param SortGenes Should the final data frame be sorted by gene name? Default is FALSE.
#' @return A new merged data frame with all the search terms combined from the lists supplied. If sort.gene=TRUE, genes will be sorted by name.
#' @description 
#' This function merges two data frames with search terms. This allows users to easily add search terms to data frames (either their
#' own or ones included in this package using data() as GenBank annotations for the same genes may vary in gene name. 
#' @examples
#' #load the list of search terms for mitochondrial genes
#' data(mtDNAterms) 
#' 
#' #Make a data.frame of new terms to add.
#' #This is a dummy example for a non-real annotation of COI, but lets pretend it is real.
#' add.name <- data.frame("COI","CDS", "CX1")
#' 
#' # make the column names the same for combination.
#' colnames(add.name) <- colnames(mtDNAterms)
#' 
#' #Run the merge search term function without sorting based on gene name.
#' new.terms <- MergeSearchTerms(add.name, mtDNAterms, SortGenes=FALSE)
#' 
#' #Run the merge search term function with sorting based on gene name.
#' new.terms <- MergeSearchTerms(add.name, mtDNAterms, SortGenes=TRUE)
#' 
#' #Merge search terms and create an additional column for introns and/or exons to extract
#' #In this example, add the trnK intron to the terms
#' 
#' ###Example With matK CDS and addint introns/exons for trnK###
#' #Subset out matK from cpDNAterms
#' cds.terms <- subset(cpDNAterms,cpDNAterms$Feature=="matK")
#' #Create a vecotr of NA so we can merge with the search terms for introns and exons
#' cds.terms <- cbind(cds.terms,(rep(NA,length(cds.terms$Feature))))
#' colnames(cds.terms)[4] <- "IntronExonNumber"
#' 
#' #Prepare a search term table for the intron and exons to remove
#' #We can start with the cpDNAterms for trnK
#' IntronExon.terms<-subset(cpDNAterms,cpDNAterms$Feature=="trnK")
#' 
#' #As we want to go for two exons, we will want the synonyms repeated as we are doing and intron
#' #and an exon
#' IntronExon.terms<-rbind(IntronExon.terms,IntronExon.terms)#duplicate the terms
#' 
#' #rep the sequence type we want to extract
#' IntronExon.terms$Type <- rep(c("intron","intron","exon","exon"))
#' IntronExon.terms$Feature <- rep(c("trnK_Intron","trnK_Exon2"),each=2)
#' IntronExon.terms <- cbind(IntronExon.terms,rep(c(1,1,2,2)))#Add intron/exon number info
#' 
#' #change column name for number info for IntronExon name
#' colnames(IntronExon.terms)[4] <- "IntronExonNumber"
#' 
#' #We can then merge everything together with MergeSearchTerms terms
#' IntronExonExampleTerms <- MergeSearchTerms(IntronExon.terms,cds.terms)
#' @export

MergeSearchTerms<-function(..., SortGenes=FALSE){
  dots <- list(...)
  check.names<-c("Feature","Type","Name","IntronExonNumber")
  if(length(unique(sapply(lapply(dots, dim), "[", 2)))>1){
    stop(paste("The input data frames are of different dimensions"))
  }
  for (i in sequence(length(dots))) {
    working <- TRUE
    if(!inherits(dots[[i]], "data.frame")) {
      stop(paste("The input object is the wrong format", sep=""))
      #working <- FALSE
    }
    if(dim(dots[[i]])[2]!=3&&dim(dots[[i]])[2]!=4) {
      stop(paste("The input object is the wrong dimensions", sep=""))
      #working <- FALSE
    }
    if(length(grep(TRUE,colnames(dots[[i]])!=check.names[1:unique(sapply(lapply(dots, dim), "[", 2))])>0)) {
      stop(paste("The columns for table",i, "do not have the proper names ('Feature','Type','Name' and if extracting introns/exons, the additional column 'IntronExonNumber')", sep=" "))
      #working <- FALSE
    }
  }
  new.terms <- rbind(...)
  if(SortGenes==TRUE) {
    new.terms <- new.terms[order(new.terms$Feature),]
  }
  row.names(new.terms)<-1:dim(new.terms)[[1]]
  return(new.terms)
}
