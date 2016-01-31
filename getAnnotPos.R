library(stringr)
library(seqinr)

accessions<-c("KJ914664","KC633221")

genes<-read.csv("gene.csv", header=TRUE)#test data
#genes<-genes[-3,]#this was just to see if it input NA

get.seq.pos<-function(accessions, genes, bank="genbank"){
choosebank(bank)#choose bank so it could be genbank or EMBL or others supported?
unique.gene.names <- unique(genes$gene)#unique gene names
seq.col.id<-paste(rep(as.vector(unique.gene.names),1,each=2),c("start","stop"),sep = ".")#these will be the column names for gene id
boundaries <- data.frame(matrix(nrow=length(accessions), ncol=2+2*length(unique.gene.names)))
#for each accession
for(i in sequence(length(accessions))){
  new.access<-strsplit(accessions[i],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
  rec<-query(paste("AC=",new.access,sep=""))#get the genbank record. Getting error related to paste and it won't show accession, but it is working.
  current.annot<-getAnnot(rec$req[[1]],nbl=2000)#I think nbl is ok, but maybe we should up it to something ridiculous just to be safe
  kill.comp<-gsub("complement\\("," ",current.annot)#kill complement()
  kill.comp<-gsub("\\)"," ",kill.comp)#kill trailing ) after complemet
  fix.Leu<-sub("tRNA-Leu|trnL","tRNA-Leu1",kill.comp, perl=TRUE)#rename leucine record.
  new.annot<-sub("tRNA-Ser|trnS","Ser1",fix.Leu)#rename serine
  spec.name<- subset(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+")),!(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+"))=="character(0)"))#get species name
  boundaries[i, 1] <- spec.name#add species name to final table
  boundaries[i, 2] <- accessions[i]#add accession number to the final table
  #for each gene
  for (j in sequence(length(unique.gene.names))) {
    genes.local <- subset(genes, genes$gene==unique.gene.names[j])#subset genes based on column 1 ID
    found.result <- matrix(nrow=1, ncol=3)#make empty result so if nothing is found, it is NA
    #for each search term combo
    for (k in sequence(dim(genes.local)[1])) {
      found.result.match <- str_match_all(paste(new.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s+/",genes.local[k,3],"=\\\"", genes.local[k,4], "\\\"", sep=""))#Match all cases for genes with duplicates tRNA in this case
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
colnames(boundaries)<-c("species","Accession",seq.col.id)#Add column names made above. Include
boundaries
}


#test why it is not working
test.string<-"Hello, I said Hello"
sub.res<-sub("Hello","qqqqqqq",test.string)#sub works down here, so I don't know
