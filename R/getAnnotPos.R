#' Finds annotation positions based on search terms that can later be used to seperate sequences into their annotated components.
#' @param accessions A vector of GenBank accession numbers.
#' @param genes A data frame of search terms. Pre-compiled search term lists are available as data with this package for mitogenomes and rDNA.
#' @param duplicate.genes A character vector containing duplicate gene names found in the annotation. Ex. serine and leucine in mitogenomes. default is NULL.
#' @return A table of start and stop positions of class annotPos for all the genes specified for all accession numbers that can be used to bust sequences using AnnotationBustR. Negative numbers indicate the sequence is a complement
#' @examples
#' ncbi.accessions<-c("FJ706343","FJ706292")#vector of two NCBI accession numbers to get the annotation positions of.
#' data(rDNAterms)#load rDNA search terms as these are records of rDNA sequences
#' my.seq.pos<-GetSeqPos(ncbi.accessions, rDNAterms, duplicate.genes= NULL)#Get rDNA gene positions for each sequence. There are no gene duplicates.
#' @export

GetSeqPos<-function(accessions, genes, duplicate.genes =NULL){
  seqinr::choosebank("genbank")#choose bank so it could be genbank or EMBL or others supported?
  unique.gene.names <- unique(genes$gene)#unique gene names
  {if (is.null(duplicate.genes))seq.col.id<-paste(rep(as.vector(unique.gene.names),1,each=2),c("start","stop"),sep = ".")#these will be the column names for gene id
    else seq.col.id<-paste(rep(c(as.vector(unique.gene.names),paste0(duplicate.genes,1)),1,each=2),c("start","stop"),sep=".")}
  boundaries <- data.frame(matrix(nrow=length(accessions), ncol=2+(2*length(unique.gene.names))+2*length(duplicate.genes)))#Make columns for the two info cats, start stop position and 
  #for each accession
  for(i in sequence(length(accessions))){
    new.access<-strsplit(accessions[i],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
    rec<-seqinr::query(paste("AC=",new.access,sep=""))#get the genbank record. Getting error related to paste and it won't show accession, but it is working.
    current.annot<-seqinr::getAnnot(rec$req[[1]],nbl=20000)#I think nbl is ok, but maybe we should up it to something ridiculous just to be safe
    genes.modified <- genes
    for (duplicate.gene.index in sequence(length(duplicate.genes))) {
      matching.terms <- unique(genes[which(genes$gene %in% duplicate.genes[duplicate.gene.index]),]$term3)
      matching.lines <- c()
      matching.genes <- genes[which(genes$gene %in% duplicate.genes[duplicate.gene.index]),]
      matching.genes$gene <- paste(matching.genes$gene, "1", sep="")
      matching.genes$term3 <- matching.genes$gene
      matching.genes <- unique(matching.genes)
      genes.modified <- rbind(genes.modified, matching.genes)
      for (term.index in sequence(length(matching.terms))) {
        matching.lines <- which(grepl(matching.terms[term.index], current.annot))
        if(length(matching.lines)>0) {
          current.annot[matching.lines[1]]<-sub(matching.terms[term.index],paste(duplicate.genes[duplicate.gene.index], "1", sep=""), current.annot[matching.lines[1]]) #sub 1st  instance
          break()
        }
      }
    }
    genes <- genes.modified
    unique.gene.names<-unique(genes$gene)
    spec.name<- subset(gsub("  ORGANISM  ", "",stringr::str_extract_all(current.annot, "  ORGANISM  \\D+")),!(gsub("  ORGANISM  ", "",stringr::str_extract_all(current.annot, "  ORGANISM  \\D+"))=="character(0)"))#get species name
    boundaries[i, 1] <- spec.name#add species name to final table
    boundaries[i, 2] <- accessions[i]#add accession number to the final table
    #for each gene
    for (j in sequence(length(unique.gene.names))) {
      genes.local <- subset(genes, genes$gene==unique.gene.names[j])#subset genes based on column 1 ID
      #found.result <- matrix(nrow=1, ncol=3)#make empty result so if nothing is found, it is NA
      #for each search term combo do regular, break if found, or try complement, followed by a break if found.
      for (k in sequence(dim(genes.local)[1])) {
        complement=FALSE
        #found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*\[\\"|\S]", sep=""))#Match all cases for genes with duplicates tRNA in this case
        found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*\\\"", sep=""))#Match all cases for genes with duplicates tRNA in this case
        {if (length(found.result.match[[1]])>0)break
          else found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+complement\\((\\d+)..(\\d+)\\)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*", sep=""))#Match all cases for genes with duplicates tRNA in this case
          {if(length(found.result.match[[1]])>0) { 
            complement=TRUE
            break}}
          #try complement
          found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+complement\\((\\d+)..(\\d+)\\)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*", sep=""))#Match all cases for genes with duplicates tRNA in this case
          if(length(found.result.match[[1]])>0) {
            complement=TRUE
            break
          }
        }
      }
      found.result <- found.result.match[[1]]#subset results
      boundaries[i, 1 + 2*j] <-ifelse(complement, as.numeric(found.result[,2])*-1, found.result[,2])#write the starting position
      boundaries[i, 2 + 2*j] <-ifelse(complement, as.numeric(found.result[,3])*-1, found.result[,3])#write the starting position
    }
    
  }
  colnames(boundaries)<-c("Species","Accession",seq.col.id)#Add column names made above. Include
  class(boundaries)<-append(class(boundaries),"Annot.Pos")#make object boundaries have class of Annot.Pos
  boundaries
}

looks<-GetSeqPos(accessions = tester, genes = mtDNAterms,duplicate.genes =c("tRNA_Ser", "tRNA_Leu"))


found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\s+""\\\"*\\\"", sep=""))#Match all cases for genes with duplicates tRNA in this case

found.result.match <- stringr::str_match_all(paste(current.annot, collapse=" "), paste(genes.local[k,2],"\\s+(\\d+)..(\\d+)\\s*+/*",genes.local[k,3],"=*\\\"*", genes.local[k,4], "\\\"*\\\"", sep=""))#Match all cases for genes with duplicates tRNA in this case
