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
#' #This is a dummy example for a non-real annoation of COI, but lets pretend it is real.
#' add.name<-data.frame("COI","CDS", "CX1")
#' 
#' # make the column names the same for combination.
#' colnames(add.name)<-colnames(mtDNAterms)
#' 
#' #Run the merge search term function without sorting based on gene name.
#' new.terms<-MergeSearchTerms(add.name, mtDNAterms, SortGenes=FALSE)
#' 
#' #Run the merge search term function with sorting based on gene name.
#' new.terms<-MergeSearchTerms(add.name, mtDNAterms, SortGenes=TRUE)
#' 
#' #Merge search terms and create an additional column for introns and/or exons to
#' #In this example, add the trnK intron to the terms
#' #create empty IntornExonNumber column for non-intron/exons
#' cp.terms<-cbind(cpDNAterms,rep(NA,length(cpDNAterms$Name)))
#' colnames(cp.terms)[4]<-"IntronExonNumber"#Name the column IntronExonNumber
#' trnK.intron.terms<-subset(cpDNAterms,cpDNAterms$Locus=="trnK")#subset trnK
#' #Create a vector of 1's the same length as the number of rows for trnK
#' trnK.terms<-cbind(trnK.intron.terms,rep(1,length(trnK.intron.terms$Name)))
#' colnames(trnK.terms)[4]<-"IntronExonNumber"#Name the column IntronExonNumber
#' #Use MergeSearchTerms to merge the modified cpDNAterms and new intron terms
#' all.terms<-MergeSearchTerms(cp.terms,trnK.terms)
#' @export

MergeSearchTerms<-function(..., SortGenes=FALSE){
  dots <- list(...)
  check.names<-c(colnames(mtDNAterms),"IntronExonNumber")
  if(length(unique(sapply(lapply(dots, dim), "[", 2)))>1){
    stop(paste("The input data frames are of different dimensions"))
  }
  for (i in sequence(length(dots))) {
    working <- TRUE
    if(class(dots[[i]])!="data.frame") {
      stop(paste("The input object is the wrong format", sep=""))
      #working <- FALSE
    }
    if(dim(dots[[i]])[2]!=3&&dim(dots[[i]])[2]!=4) {
      stop(paste("The input object is the wrong dimensions", sep=""))
      #working <- FALSE
    }
    if(length(grep(TRUE,colnames(dots[[i]])!=check.names[1:unique(sapply(lapply(dots, dim), "[", 2))])>0)) {
      stop(paste("The columns for table",i, "do not have the proper names ('Locus','Type','Name' and if extracting introns/exons, the additional column 'IntronExonNumber')", sep=" "))
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
