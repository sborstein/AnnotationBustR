#' Allows for the merging of two tables containing search terms to expand search term database for getAnnotPos.
#' @param ... the data frames of search terms you want to combine into a single list. The Data frames should have stringsAsFactors=FALSE if you want to sort them.
#' @param sort.gene Should the final data frame be sorted by gene name? Default is FALSE.
#' @return A new merged data frame with all the search terms combined from the lists supplied. If sort.gene=TRUE, genes will be sorted by name.
#' @export

MergeSearchTerms<-function(..., sort.gene=FALSE){
  dots <- list(...)
  for (i in sequence(length(dots))) {
    working <- TRUE
    if(dim(dots[[i]])[2]!=4) { 
      working <- FALSE
    }
    if(class(dots[[i]])!="data.frame") {
      working <- FALSE
    }
    if (!working) {
      stop(paste("The ", i, "th input object is the wrong format", sep=""))
    }
  }
  new.terms <- rbind(...)
  if(sort.gene) {
    new.terms <- new.terms[order(new.terms$gene),]
  }
 
  return(new.terms)
}
