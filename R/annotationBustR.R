#' Breaks up genbank sequences into their annotated components based on positions found using the get.seq.pos function.
#' @param SeqPos Object of class Annot.Pos with the starting and stopping position for all genes.
#' @return Writes a fasta file(s) to the working directory for each gene in the object of class Annot.Pos supplied into the function.
#' @examples
#' ncbi.accessions<-c("FJ706343","FJ706292")#vector of two NCBI accession numbers to get the annotation positions of.
#' data(rDNAterms)#load rDNA search terms as these are records of rDNA sequences
#' my.seq.pos<-GetSeqPos(ncbi.accessions, rDNAterms, duplicate.genes= NULL)#Get rDNA gene positions for each sequence. There are no gene duplicates.
#' AnnotationBust(my.seq.pos)#Run the function to extract genes using positions found in GetSeqPos and write them to FASTA file in the working directory
#' @export

AnnotationBust<-function(SeqPos){
if(length((class(SeqPos))) <2)#check to make sure class is length 2
  stop("Input is not a data.frame or of class Annot.Pos")
if(!(class(SeqPos)[2]=="Annot.Pos"))#if it is 2, check that it is Annot.Pos
  stop("Input is not a data.frame or of class Annot.Pos")
gene.starts <- which(grepl("start", colnames(SeqPos)))#get the start position for each sequence
for (gene.index in sequence(length(gene.starts))) {
  local.info <- data.frame(Species=SeqPos$Species, Accession=SeqPos$Accession, Start=SeqPos[,gene.starts[gene.index]], Stop=SeqPos[,1+gene.starts[gene.index]], stringsAsFactors=FALSE)#For each start, pul out the relevant info and next column which is the stop. Make it a non-factor.
  local.info<- subset(local.info, is.na(local.info[,3])==FALSE)#For the current gene, remove entries from the column that are NA.
  gene.name<-strsplit(colnames(SeqPos)[gene.starts[gene.index]],".",fixed=TRUE)[[1]][1]
  open.command="w"
  for (taxon.index in sequence(dim(local.info)[1])) {
    current.seq<-ape::read.GenBank(local.info$Accession[taxon.index],species.names=TRUE,as.character=TRUE)#access number for read genbank
    if (local.info[taxon.index,3]<1){#check if it is a complement, as it will start with a "-"
      cut.seq<-current.seq[[1]][(local.info$Start[taxon.index]*-1):(local.info$Stop[taxon.index]*-1)]#if negative multiply by one to cut based on start/stop position
      cut.seq<-rev(comp(cut.seq))}#make the reverse complement as it is a complementary seq
     else {cut.seq<-current.seq[[1]][local.info$Start[taxon.index]:local.info$Stop[taxon.index]]}#cut based on start/stop position
    #cut.seq<-current.seq[[1]][local.info$Start[taxon.index]:local.info$Stop[taxon.index]]#cut based on start/stop position
    seq.ids<-attr(current.seq, "species")#extract attribute-species name, from genbank takedown above
    seqinr::write.fasta(sequences=cut.seq, names=seq.ids, file.out=paste(gene.name,"fa", sep="."),open=open.command)#write to fasta with gene name replacing the accession# with the species name
    open.command="a"
    }
  }
}
