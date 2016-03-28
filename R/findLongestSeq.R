#' Find the longest sequence for each species from a list of GenBank accession numbers.
#' 
#' @param accessions A vector of GenBank accession numbers.
#' @return A list of genbank accessions numbers for the longest sequence for each taxon in a list of accession numbers.
#' @example 
#' genbank.accessions<-c("KP978059.1","KP978060.1","JX516105.1","JX516111.1")#vector of 4 genbank accessions, two each for two species
#' long.seq.result<-Find.Longest.Seq(genbank.accessions)#returns the longest sequence respectively for the two species.

Find.Longest.Seq<-function(accessions){
  multi.ncbi.hits<-NULL#empty vector for merging records
  final.accession<-NULL#empty vector
  x<-seq_along(numbs)#seq along all to get number
  split.access <- split(numbs, ceiling(x/600))#break them up into groups of 600
  for (i in 1:length(split.access)){
    ncbi.hits<-ncbi_byid(split.access[[i]])[,c(1,3:5)]#use the ropensci traits package to read the accession numbers
    multi.ncbi.hits<-rbind(multi.ncbi.hits,ncbi.hits)#merge together all search hits given by merging to the ith multiple of 600
  }
  unique.taxa<-unique(multi.ncbi.hits$taxon)#get the names of taxa to find the longest seq
  for(j in sequence(length(unique.taxa))){#for each taxa
    current.spec<-subset(multi.ncbi.hits,multi.ncbi.hits$taxon==unique.taxa[[j]])#get the data for all the records for species j, minus the sequence
    longest.seq<-subset(current.spec,current.spec$length==sort(as.numeric(current.spec$length),decreasing = TRUE)[1])[1,]#for all of the
    final.accession<-rbind(final.accession,longest.seq)#combine individual rows into a data frame
  }
  row.names(final.accession)<-1:dim(final.accession)[1]#make row names the length of 
  final.accession#return final list
}