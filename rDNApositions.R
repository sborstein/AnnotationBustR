#' Gets the location of individual genes in a rDNA  Sequence and produces an object of class rDNA.locs
#' @param accessions A vector of GenBank accession numbers. Must not contain RefSeq number prefix, (i.e. NC_,XM_, etc.)
#' @param bank Database to access accession number. Default is genbank. Options follow those of seqinr. 

rDNA.Locs<-function(accessions, bank="genbank"){
  choosebank(bank)#choose bank so it could be genbank or EMBL or others supported
  #empty vectors for filling with position info
  rRNA18s<-c(NULL,NULL)
  its.1<-c(NULL,NULL)
  rRNA5.8s<-c(NULL,NULL)
  its.2<-c(NULL,NULL)
  rRNA28s<-c(NULL,NULL)
  organism<-c(NULL,NULL)
  for(i in sequence(length(accessions))){
    new.access<-strsplit(accessions[i],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
    rec<-query(paste0("AC=",new.access))#paste in accession number for i
    current.annot<-getAnnot(rec$req[[1]],nbl=2000)
    kill.comp<-gsub("<","",current.annot)
    new.annot<-gsub(">","",kill.comp)
    spec.names<- subset(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+")),!(gsub("  ORGANISM  ", "",str_extract_all(new.annot, "  ORGANISM  \\D+"))=="character(0)"))#get species name
    rRNA.18s<- str_match(paste(new.annot, collapse=" "), paste("rRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "18S ribosomal RNA", "\\\"", sep=""))
    ITS1<-str_match(paste(new.annot, collapse=" "), paste("misc_RNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "internal transcribed spacer 1", "\\\"", sep=""))
    rRNA.5.8s<- str_match(paste(new.annot, collapse=" "), paste("rRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "5.8S ribosomal RNA", "\\\"", sep=""))
    ITS2<-str_match(paste(new.annot, collapse=" "), paste("misc_RNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "internal transcribed spacer 2", "\\\"", sep=""))
    rRNA.28s<- str_match(paste(new.annot, collapse=" "), paste("rRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "28S ribosomal RNA", "\\\"", sep=""))
    #Append lists
    organism<-rbind(organism,spec.names[1])
    accession<-accessions
    rRNA18s<-rbind(rRNA18s,as.numeric(rRNA.18s[,2:3])) 
    its.1<-rbind(its.1,as.numeric(ITS1[,2:3]))
    rRNA5.8s<-rbind(rRNA5.8s,as.numeric(rRNA.5.8s[,2:3]))
    its.2<-rbind(its.2,as.numeric(ITS2[,2:3]))
    rRNA28s<-rbind(rRNA28s,as.numeric(rRNA.28s[,2:3]))
  }
  all.rDNA.locs<-cbind(rRNA18s,its.1,rRNA5.8s,its.2,rRNA28s)
  colnames(all.rDNA.locs)<-c("rRNA18s.start","rRNA18s.stop","ITS1.start","ITS1.stop","rRNA5.8s.start","rRNA5.8s.stop","ITS2.start","ITS2.stop","rRNA28s.start","rRNA28s.stop")
  all.rDNA.locs.info<-cbind(accessions,all.rDNA.locs)
  row.names(all.rDNA.locs.info)<-organism
  all.rDNA.locs.info
}


###test####
gb.numbers<-c("AB277057.1","FJ706343.1")#one of these is incomplete, should fill in with NA as it won't get hits
my.pos<-rDNA.Locs(gb.numbers, bank="genbank")#Works
