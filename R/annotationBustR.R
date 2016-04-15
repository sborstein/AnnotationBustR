library(ape)
library(seqinr)

#' Breaks up genbank sequences into their annotated components based on positions found using the get.seq.pos function.
#' @param Seq.Positions Object of class Annot.Pos with the starting and stopping position for all genes.
#' @return Returns a fasta file for each gene in the object of class Annot.Pos supplied into the function.
#' @export

AnnotationBust<-function(Seq.Positions){
gene.starts <- which(grepl("start", colnames(Seq.Positions)))#get the start position for each sequence
for (gene.index in sequence(length(gene.starts))) {
  local.info <- data.frame(Species=Seq.Positions$Species, Accession=Seq.Positions$Accession, Start=Seq.Positions[,gene.starts[gene.index]], Stop=Seq.Positions[,1+gene.starts[gene.index]], stringsAsFactors=FALSE)#For each start, pul out the relevant info and next column which is the stop. Make it a non-factor.
  local.info <- local.info[-which(is.na(local.info[,3])),]#For the current gene, remove entries from the column that are NA.
  #gene.name <- gsub('.start', colnames(Seq.Positions)[gene.starts[gene.index]])
  gene.name<-strsplit(colnames(Seq.Positions)[gene.starts[gene.index]],".",fixed=TRUE)[[1]][1]
  open.command="w"
  for (taxon.index in sequence(dim(local.info)[1])) {
    current.seq<-ape::read.GenBank(local.info$Accession[taxon.index],species.names=TRUE,as.character=TRUE)#access number for read genbank
    cut.seq<-current.seq[[1]][local.info$Start[taxon.index]:local.info$Stop[taxon.index]]#cut based on start/stop position
    seq.ids<-attr(current.seq, "species")#extract attribute-species name, from genbank takedown above
    seqinr::write.fasta(sequences=cut.seq, names=seq.ids, file.out=paste(gene.name,"fa", sep="."),open=open.command)#write to fasta with gene name replacing the accession# with the species name
    open.command="a"
    }
  }
}