#' Converts object of class Annot.Pos into an accession table
#' @param SeqPos Object that is a data frame of class Annot.Pos returned by the function GetSeqPos
#' @return Returns a data frame with rowss as species names and columns as loci filled with the accessions used for each species/gene.
#' @examples
#' ncbi.accessions<-c("FJ706343","FJ706292")#vector of two NCBI accession numbers to get the annotation positions of.
#' data(rDNA.Genes)#load rDNA search terms
#' my.seq.pos<-GetSeqPos(ncbi.accessions, rDNA.Genes, duplicate.genes= NULL)#Get rDNA gene positions for each sequence. There are no gene duplicates.
#' my.accessions<-MakeAccessionTable(SeqPos=my.seq.pos)#Convert my.seq.pos, which is the class data.frame and Annot.Pos into a table of accessions for each gene.
#'  @export

MakeAccessionTable<-function(SeqPos){
  if(length((class(SeqPos))) <2)#check to make sure class is length 2
    stop("Input is not a data.frame or of class Annot.Pos")
  if(!(class(SeqPos)[2]=="Annot.Pos"))#if it is 2, check that it is Annot.Pos
    stop("Input is not a data.frame or of class Annot.Pos")
  UniqueSpecies<-unique(SeqPos$Species)
  gene.starts <- which(grepl("start", colnames(SeqPos)))#get the start position for each sequence
  local.info <- data.frame(Species=SeqPos$Species, Accession=SeqPos$Accession, SeqPos[,gene.starts], stringsAsFactors=FALSE)#For each start, pul out the relevant info and next column which is the stop. Make it a non-factor.
  present<-!is.na(local.info[,3:dim(local.info)[2]])#is.NA is FALSE if filled with number
  local.info[3:dim(local.info)[2]]<-present#Merge the true/false
  add.access<-ifelse(local.info[,3:dim(local.info)[2]]==TRUE, local.info$Accession, NA)#Paste in accessions
  local.info[3:dim(local.info)[2]]<-add.access#merge in the list of accessions
  local.info<-local.info[,-2]#kill accession number as it is already filled in the table
  gene.names<-gsub(".start", "",colnames(local.info))
  accession.table<-data.frame(matrix(nrow=length(UniqueSpecies), ncol=length(local.info)))#make accession table length of species * length of loci
  colnames(accession.table)<-gene.names
  #For each species get the accessions for each gene
  for (i in 1:length(UniqueSpecies)){
    current.spec<-subset(local.info,local.info$Species==UniqueSpecies[i])
    accession.table[i,1]<-UniqueSpecies[i]
    #start at 2 for column loop as first column is species name and put into accession table using the line above
    #For each gene, get the numbers if they are present
    for (j in 2:length(local.info)){
      current.loci<-current.spec[,j]#The current column/loci
      accessions<-subset(current.loci,!is.na(current.loci))#subset the ones that are present
      numbers<-ifelse(length(accessions)==0, NA, paste(accessions, sep=",", collapse=","))#ifelse statemen, if not present, fill with NA
      accession.table[i,j]<-numbers
      }
    }
    accession.table
}
