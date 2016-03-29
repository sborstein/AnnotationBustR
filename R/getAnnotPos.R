library(stringr)
library(seqinr)

#' Finds annotation positions based on search terms that can later be used to seperate sequences into their annotated components.
#' @param accessions A vector of GenBank accession numbers.
#' @param genes A data frame of search terms. Pre-compiled search term lists are available as data with this package for mitogenomes and rDNA.
#' @param bank Name of bank, either genbank or embl. Default is genbank.
#' @return A table of start and stop positions of class annotPos for all the genes specified for all accession numbers that can be used to bust sequences using AnnotationBustR.

get.seq.pos<-function(accessions, genes, bank="genbank"){
choosebank(bank)#choose bank so it could be genbank or EMBL or others supported?
unique.gene.names <- unique(genes$gene)#unique gene names
seq.col.id<-paste(rep(as.vector(unique.gene.names),1,each=2),c("start","stop"),sep = ".")#these will be the column names for gene id
boundaries <- data.frame(matrix(nrow=length(accessions), ncol=2+2*length(unique.gene.names)))
#for each accession
for(i in sequence(length(accessions))){
  new.access<-strsplit(accessions[i],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
  rec<-query(paste("AC=",new.access,sep=""))#get the genbank record. Getting error related to paste and it won't show accession, but it is working.
  current.annot<-getAnnot(rec$req[[1]],nbl=20000)#I think nbl is ok, but maybe we should up it to something ridiculous just to be safe
  new.annot<-gsub("complement\\(|\\(|<"," ",current.annot)#kill complement()
  new.annot<-gsub("\\)|>","",new.annot)#kill trailing ) after complemet
  matching.lines.Leu<-which(grepl("tRNA-Leu|trnL|trnL-uaa|trnL TAA",new.annot))#find Leucine lines
  new.annot[matching.lines.Leu[1]]<-sub("tRNA-Leu|trnL|trnL-uaa|trnL TAA","tRNA-Leu1",new.annot[matching.lines.Leu[1]])#sub 1st leucine instance
  matching.lines.Ser<-which(grepl("tRNA-Ser|trnS|trnS-uga|trnS TGA",new.annot))#find serine lines
  new.annot[matching.lines.Ser[1]]<-sub("tRNA-Ser|trnS|trnS-uga|trnS TGA","tRNA-Ser1",new.annot[matching.lines.Ser[1]])#sub first instance of serine
  spec.name<- subset(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+")),!(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+"))=="character(0)"))#get species name
  boundaries[i, 1] <- spec.name#add species name to final table
  boundaries[i, 2] <- accessions[i]#add accession number to the final table
  #for each gene
  for (j in sequence(length(unique.gene.names))) {
    genes.local <- subset(genes, genes$gene==unique.gene.names[j])#subset genes based on column 1 ID
    found.result <- matrix(nrow=1, ncol=3)#make empty result so if nothing is found, it is NA
    #for each search term combo
    for (k in sequence(dim(genes.local)[1])) {
      found.result.match <- str_match_all(paste(new.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*", sep=""))#Match all cases for genes with duplicates tRNA in this case
      #if statement to break searching when a result is found
      if((dim(found.result.match[[1]])[1])>0) {
        found.result <- found.result.match[[1]][1,]#subset results
        boundaries[i, 1 + 2*j] <- found.result[2]#write the starting position
        boundaries[i, 2 + 2*j] <- found.result[3]#write stop position
        break
      }
    }
  }
}
colnames(boundaries)<-c("Species","Accession",seq.col.id)#Add column names made above. Include
class(boundaries)<-append(class(boundaries),"Annot.Pos")#make object boundaries have class of Annot.Pos
boundaries
}
