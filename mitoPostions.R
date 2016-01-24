#' Gets the location of individual genes in a Mitochondrial Genome Sequence
#' @param accessions A vector of GenBank accession numbers
#' @param bank database to access accession number. Default is genbank. Options follow those of seqinr.

Mito.Locs<-function(accessions, bank="genbank"){
  choosebank(bank)#choose bank so it could be genbank or EMBL
  for(i in sequence(length(accessions))){
    seq<-accessions[[i]]#This chokes hard for some reason and makes the annotations funky
    current.record<-query(paste("AC=",seq, sep=""))
    current.annot<-getAnnot(current.record$req,nbl=2000)

    
  }
}


current.annot
new.annot<-gsub("complement\\("," ",current.annot)
