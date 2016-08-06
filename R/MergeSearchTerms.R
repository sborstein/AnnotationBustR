#' Allows for the merging of two tables containing search terms to expand search term database for getAnnotPos.
#' @param ... the data frames of search terms you want to combine into a single data frame The Data frame(s) should have stringsAsFactors=FALSE listed if you want to sort them.
#' @param sort.gene Should the final data frame be sorted by gene name? Default is FALSE.
#' @return A new merged data frame with all the search terms combined from the lists supplied. If sort.gene=TRUE, genes will be sorted by name.
#' @description 
#' This function merges two data frames with search terms. This allows users to easily add search terms to data frames (either their
#' own or ones included in this package using data()) as GenBank annotations for the same gene may vary in gene name. 
#' @examples
#' data(mtDNAterms) #load the list of search terms for mitochondrial genes
#' add.name<-data.frame("COI","CDS", "CX1")#Make a data.frame of new terms to add. This is a dummy example for a non-real annoation of COI, but lets pretend it is real.
#' colnames(add.name)<-colnames(mtDNAterms)# make the column names the same for combination.
#' MergeSearchTerms(add.name, mtDNAterms, SortGenes=FALSE)#Run the merge search term function without sorting based on gene name.
#' MergeSearchTerms(add.name, mtDNAterms, SortGenes=TRUE)#Run the merge search term function with sorting based on gene name.
#' @export

MergeSearchTerms<-function(..., SortGenes=FALSE){
  dots <- list(...)
  for (i in sequence(length(dots))) {
    working <- TRUE
    if(class(dots[[i]])!="data.frame") {
      stop(paste("The input object is the wrong format", sep=""))
      #working <- FALSE
    }
    if(dim(dots[[i]])[2]!=3) {
      stop(paste("The input object is the wrong dimensions", sep=""))
      #working <- FALSE
    }
  }
  new.terms <- rbind(...)
  if(SortGenes==TRUE) {
    new.terms <- new.terms[order(new.terms$Locus),]
  }
  row.names(new.terms)<-1:dim(new.terms)[[1]]
  return(new.terms)
}
