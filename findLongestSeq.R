library(traits)
library(ape)

#' Find the longest sequence for each species from a list of GenBank accession numbers.
#' @param accessions A vector of GenBank accession numbers.
#' @return A list of genbank accessions numbers for the longest sequence for each taxon in a list of accession numbers.

Find.Longest.Seq<-function(accessions){
      ncbi.hits<-ncbi_byid(accessions, format="fasta")[,c(1,3:5)]#use the ropensci traits package to read the accession numbers
      get.names<-read.GenBank(accessions,species.name=TRUE)#use ape read genbank for getting around trinomial issue
      Species = attr(get.names, "species")#get species from ape
      ncbi.hits$taxon<-Species#sub the ape names in for the 
      unique.taxa<-unique(ncbi.hits$taxon)#get the names of taxa to find the longest seq
      final.accession<-NULL
      for(i in sequence(length(unique.taxa))){#for each taxa
        current.spec<-subset(ncbi.hits,ncbi.hits$taxon==unique.taxa[[i]])#get the data for all the records for species i, minus the sequence
        longest.seq<-subset(current.spec,current.spec$length==sort(current.spec$length,decreasing = TRUE)[1])[1,]#for all of the
        final.accession<-rbind(final.accession,longest.seq)#combine individual rows into a data frame
      }
      row.names(final.accession)<-1:dim(final.accession)[1]#make row names the length of 
      final.accession#return final list
    }