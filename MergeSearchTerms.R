#' Allows for the merging of two tables containing search terms to expand search term database for getAnnotPos.
#' @param search.lists A list containing the data frames of search terms you want to combine into a single list.
#' @param sort.loci Should the final data frame be sorted by gene name? Default is FALSE.
#' @return A new merged data frame with all the search terms combined from the lists supplied.

list1<-read.csv("MitoGenesList.csv")
list2<-read.csv("rDNA.csv")
test.list<-list(list1,list2)

MergeSearchTerms<-function(search.lists, sort.gene=FALSE){
  new.list<-ldply(search.lists, data.frame)
  if(sort.gene==TRUE){
    new.list<-new.list[order(new.list$gene),]#not fully working, sorts, but the original data frames are still on top of each other.
    new.list
  }
  else{
    new.list
  }
}
