#' Gets the location of individual genes in a Mitochondrial Genome Sequence
#' @param accessions A vector of GenBank accession numbers. Must not contain RefSeq number prefix, (i.e. NC_)
#' @param bank Database to access accession number. Default is genbank. Options follow those of seqinr. 

Mito.Locs<-function(accessions, bank="genbank"){
  choosebank(bank)#choose bank so it could be genbank or EMBL or others supported?
  Phe<-c(NULL,NULL)
  for(i in sequence(length(accessions))){
    new.access<-strsplit(accessions[i],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
    rec<-query(paste("AC=",new.access,sep=""))
    current.annot<-getAnnot(rec$req[[1]],nbl=2000)
    kill.comp<-gsub("complement\\("," ",current.annot)
    new.annot<-gsub("\\)"," ",kill.comp)
    #start pulling out parts of mitogenome, should be pretty obvious what each part does
    tRNA.Phe<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Phe", "\\\"", sep=""))
    rRNA.12s<- str_match(paste(new.annot, collapse=" "), paste("rRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "12S ribosomal RNA", "\\\"", sep=""))
    tRNA.Val<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Val", "\\\"", sep=""))
    rRNA.16s<- str_match(paste(new.annot, collapse=" "), paste("rRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "16S ribosomal RNA", "\\\"", sep=""))
    tRNA.Leu<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Leu", "\\\"", sep=""))
    tRNA.Leu1<-tRNA.Leu[[1]][1,]#subset first instance of tRNA-Leu gene
    tRNA.Leu1<-tRNA.Leu1[2:3]#subset first instance of tRNA-Leu gene
    tRNA.Leu2<-tRNA.Leu[[1]][2,]#subset second instance of tRNA-Leu
    tRNA.Leu2<-tRNA.Leu2[2:3]#subset second instance of tRNA-Leu
    ND1 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND1", "\\\"", sep=""))
    tRNA.Ile<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Ile", "\\\"", sep=""))
    tRNA.Gln<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Gln", "\\\"", sep=""))
    tRNA.Met<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Met", "\\\"", sep=""))
    ND2 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND2", "\\\"", sep=""))
    tRNA.Trp<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Trp", "\\\"", sep=""))
    tRNA.Ala<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Ala", "\\\"", sep=""))
    tRNA.Asn<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Asn", "\\\"", sep=""))
    tRNA.Cys<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Cys", "\\\"", sep=""))
    tRNA.Tyr<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Tyr", "\\\"", sep=""))
    COI <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "COX1", "\\\"", sep=""))
    tRNA.Ser<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Ser", "\\\"", sep=""))
    tRNA.Ser1<-tRNA.Ser[[1]][1,]#subset first instance of tRNA-Ser gene
    tRNA.Ser1<-tRNA.Ser1[2:3]#subset first instance of tRNA-Ser gene
    tRNA.Ser2<-tRNA.Ser[[1]][2,]#subset second instance of tRNA-Ser
    tRNA.Ser2<-tRNA.Ser2[2:3]#subset second instance of tRNA-Ser
    tRNA.Asp<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Asp", "\\\"", sep=""))
    COII <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "COX2", "\\\"", sep=""))
    tRNA.Lys<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Lys", "\\\"", sep=""))
    ATP8 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ATP8", "\\\"", sep=""))
    ATP6 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ATP6", "\\\"", sep=""))
    COIII <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "COX3", "\\\"", sep=""))
    tRNA.Gly<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Gly", "\\\"", sep=""))
    ND3 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND3", "\\\"", sep=""))
    tRNA.Arg<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Arg", "\\\"", sep=""))
    ND4L <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND4L", "\\\"", sep=""))
    ND4 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND4", "\\\"", sep=""))
    tRNA.His<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-His", "\\\"", sep=""))
    ND5 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND5", "\\\"", sep=""))
    ND6 <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "ND6", "\\\"", sep=""))
    tRNA.Glu<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Glu", "\\\"", sep=""))
    CYTB <- str_match(paste(new.annot, collapse=" "), paste("gene\\s+(\\d+)..(\\d+)\\s+/gene=\\\"", "CYTB", "\\\"", sep=""))
    tRNA.Thr<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Thr", "\\\"", sep=""))
    tRNA.Pro<-str_match_all(paste(new.annot, collapse=" "), paste("tRNA\\s+(\\d+)..(\\d+)\\s+/product=\\\"", "tRNA-Pro", "\\\"", sep=""))
    Dloop <- str_match(paste(new.annot, collapse=" "), paste("D-loop", "\\s+(\\d+)..(\\d+)", sep=""))
  }
}


#Best way to capture output?
#Make Class Mito.Locs for MitoBustR