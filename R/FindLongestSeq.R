#' Find the longest sequence for each species from a list of GenBank accession numbers.
#' @param Accessions A vector of GenBank accession numbers.
#' @details For a set of GenBank accession numbers, this will return the longest sequence for in the set for species.
#' @return A list of genbank accessions numbers for the longest sequence for each taxon in a list of accession numbers.
#' @examples 
#' #a vector of 4 genbank accessions, there are two for each species.
#' genbank.accessions<-c("KP978059.1","KP978060.1","JX516105.1","JX516111.1")
#' \dontrun{
#' #returns the longest sequence respectively for the two species.
#' long.seq.result<-FindLongestSeq(genbank.accessions)
#' }
#' @export

FindLongestSeq<-function(Accessions){
Accessions<-sapply(strsplit(as.character(Accessions),"\\.",perl = TRUE), `[`, 1)
raw.accessions<-ape::read.GenBank(Accessions)#get the seqs
seq.ids<-attr(raw.accessions, "species")#extract attribute-species name, from genbank takedown above
seq.data<-data.frame(cbind(attr(raw.accessions, "species"), names(raw.accessions),as.vector(summary(raw.accessions)[,1])),stringsAsFactors = FALSE)#data frame it
colnames(seq.data)<-c("Species","Accession","Length")#attribute column names
seq.data$Length<-as.numeric(seq.data$Length)#make the length of the sequence numeric
uni.taxa<-unique(seq.data$Species)#get unique taxa
long.seqs<-data.frame(matrix(nrow=length(uni.taxa), ncol=dim(seq.data)[2]))#empty data frame to store results
colnames(long.seqs)<-colnames(seq.data)
  for (taxa.index in 1:length(uni.taxa)){
    current.tax<-subset(seq.data,seq.data$Species==uni.taxa[taxa.index])
    longest.seq<-subset(current.tax,current.tax$Length==sort(as.numeric(current.tax$Length),decreasing = TRUE)[1])[1,]#id and grab longest seq
    long.seqs[taxa.index, ]<-longest.seq#append data frame with longest seq
  }
long.seqs#return longest seqs
}
