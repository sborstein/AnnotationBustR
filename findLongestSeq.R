#' Gets sequence lengths for a list of GenBank accession numbers and returns the GenBank accesion number with the longest sequence per species
#' @param accessions A vector of GenBank accession numbers


Find.Longest.Seq<-function(accessions){
      ncbi.hits<-ncbi_byid(accessions, format="fasta")[,c(1,3:5)]#use the ropensci traits package to read the accession numbers
      unique.taxa<-unique(ncbi.hits$taxon)#get the names of taxa to find the longest seq
      final.accession<-NULL
      for(i in sequence(length(unique.taxa))){#for each taxa
        current.spec<-subset(ncbi.hits,ncbi.hits$taxon==unique.taxa[[i]])#get the data for all the records for species i, minus the sequence
        longest.seq<-subset(current.spec,current.spec$length==sort(current.spec$length,decreasing = TRUE)[1])[1,]#for all of the
        final.accession<-rbind(final.accession,longest.seq)#combine individual rows into a data frame
      }
    }
