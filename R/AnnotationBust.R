#' Breaks up genbank sequences into their annotated components based on a set of search terms and writes each sub-sequence of interest to a file for each accession number supplied.
#' @param Accessions A vector of accession numbers. Note that refseq numbers (i.e. prefixes like XM_ and NC_) will not work.
#' @param Terms A data frame of search terms. Search terms for animal mitogenomes, nuclear rRNA, chloroplast genomes, and plant mitogenomes are pre-made and can be loaded using the data()function. Additional terms can be addded using the MergeSearchTerms function.
#' @param Duplicates A vector of loci names that have more than one copy. Default is NULL
#' @param DuplicateInstances A numeric vector the length of Duplicates of the number of duplicates for each duplicated gene you wish to extract.Default is NULL
#' @param TranslateSeqs Logical as to whether or not the sequence should be translated to the corresponding peptide sequence. Default is FALSE. Note that this only makes sense to list as TRUE for protein coding sequences.
#' @param TranslateCode Numerical representing the GenBank translation code for which sequences should be translated under. Default is 1. For all options see ?seqinr::getTrans. Some possible useful ones are: 2. Vertebrate Mitochondrial, 5. Invertebrate Mitochondrial, and 11. bacterial+plantplastid
#' @param DuplicateSpecies Logical. As to whether there are duplicate individuals per species. If TRUE, adds the accession number to the fasta file 
#' @details The AnnotationBust function takes a vector of accession numbers and a data frame of search terms and extracts sub-sequences from genomes or concatenated sequences.
#' This function requires internet access. It writes files in the FASTA format to the working directory and returns an accession table. AnnoitationBustR comes with pre-made
#' search terms for mitogenomes, chloroplast genomes, and rDNA that can be loaded using data(mtDNAterms),data(cpDNAterms), and data(rDNAterms) respectively.
#' Search terms can be completely made by the user as long as they follow a similar format, with three columns. The first, Locus should contain the name of the files to be written. We recommend following
#' a similar naming convention to what we currently have in the pre-made data frames to ensure that files are named properly, characters like "-" or "." should be avoided as to not throw off R.
#' The second column, Type contains the type of sub-sequence it is, with options being CDS, rRNA, tRNA, misc_RNA, and D_Loop. The last column, Name, consists of a
#' name for the locus of interest. For numerous synonyms for the same locus, one should have each synonym as its own row. For a more detailed walkthrough on using
#' AnnotationBust you can call the vignette with vignette("AnnotationBustR).
#' @return Writes a fasta file(s) to the current working directory selected for each unique sub-sequence of interest in Terms containing all the accession numbers the sub-sequence was fond in
#' @return Writes an data.frame of the accession numbers per loci that can be turned into an accession table using the function MakeAccessionTable
#' @examples
#' ncbi.accessions<-c("FJ706295","FJ706343","FJ706292")#vector of two NCBI accession numbers to get the annotation positions of.
#' data(rDNA.Genes)#load rDNA search terms from AnnotationBustR
#' my.sequences<-AnnotationBust(ncbi.accessions, rDNA.Genes, DuplicateSpecies=TRUE)#Run AnnotationBustR and write files to working directory
#' my.sequences#Return the accession table for each species.
#' @export

AnnotationBust<-function(Accessions, Terms, Duplicates= NULL,DuplicateInstances=NULL, TranslateSeqs=FALSE, TranslateCode=1, DuplicateSpecies=FALSE){
  seqinr::choosebank("genbank")
  uni.locus<-unique(Terms$Locus)
  uni.type<-unique(Terms$Type)
  ##Deal with duplicates in regards to writing output files##
  if(is.null(Duplicates)==FALSE){
    new.file.names<-character(length = 0)
    dup.frame<-data.frame(Duplicates, DuplicateInstances)
    singles<-uni.locus[uni.locus %in% dup.frame$Duplicates==FALSE]#get matching duplicates
    doubles<-uni.locus[uni.locus %in% dup.frame$Duplicates==TRUE]#get matching duplicates
    for (i in 1:length (doubles)){
      sub.dups<-subset(dup.frame, dup.frame$Duplicates %in% doubles[i]==TRUE)#subset the duplicate genes
      number.names<-paste0(sub.dups$Duplicates, 1:sub.dups$DuplicateInstances)#add duplicate numbers to names
      new.file.names<-append(new.file.names, number.names)#append the new pasted file names
    }
    file.names<-c(as.vector(singles), new.file.names)#combo the names of single and dup loci
  }else {file.names<-as.vector(uni.locus)}
  Accession.Table<-data.frame(data.frame(matrix(NA, nrow = length(Accessions), ncol = 1+length(file.names))))#make empty accession table
  colnames(Accession.Table)<-c("Species",file.names)#attach loci names as colnames
  for (locus.index in 1:length(file.names)){#For each locus, write empty file to append to later
    write(NULL, file=paste0(file.names[locus.index], ".fasta"))#write it as a fasta
  }
  for (subsequence.type.index in 1:length(uni.type)){#for each sequence type, get the stuff subset. Don't do d-loop as it has a seperate pathway to being captured
    if (uni.type[subsequence.type.index]=="CDS"){
      CDS.Search<-subset(Terms, Terms$Type=="CDS")
      unique.CDS<-unique(CDS.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="tRNA"){
      tRNA.Search<-subset(Terms, Terms$Type=="tRNA")
      unique.tRNA<-unique(tRNA.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="rRNA"){
      rRNA.Search<-subset(Terms, Terms$Type=="rRNA")
      unique.rRNA<-unique(rRNA.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="misc_RNA"){
      misc.RNA.search<-subset(Terms, Terms$Type=="misc_RNA")
      unique.misc_RNA<-unique(misc.RNA.search$Locus)
    }
  }
  for (accession.index in 1:length(Accessions)){
    new.access<-strsplit(Accessions[accession.index],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
    species.name<-attr(ape::read.GenBank(Accessions[accession.index]),"species")#get the sequence names
    ifelse(DuplicateSpecies==TRUE, seq.name<-paste(species.name,new.access,sep = "_"), seq.name<-species.name)
    Accession.Table$Species[accession.index]<-species.name
    full.rec<-seqinr::query(paste0("AC=",new.access))
    full.annot<-seqinr::getAnnot(full.rec$req, nbl=20000)#read in the annotation
    start.loci<-grep("FEATURES", full.annot[[1]])#find start of features
    new.ann<-full.annot[[1]][-c(1:start.loci,length(full.annot[[1]]))]#cut just to the feature
    new.ann<-gsub("D-loop","Dloop",new.ann)#kill wild character "-"
    #Start parsing into a list to speed up. Only have to ping ACNUC one time this way
    annotation.list <- list()
    to.store <- c()
    for (i in sequence(length(new.ann))) {
      if(grepl("     [a-zA-Z_]\\w*.    ", new.ann[i]) & i!=1) {
        annotation.list[[length(annotation.list)+1]] <- to.store
        to.store <- new.ann[i]
      } else {
        to.store <- append(to.store, new.ann[i])
      }
    }
    annotation.list[[length(annotation.list)+1]] <- to.store
    #Now find loci for every loci type
    for (loci.type.index in 1:length(uni.type)){
      if (uni.type[loci.type.index]=="tRNA")  {
        rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=tRNA", sep=" "))#get the tRNA
        current.annot<-annotation.list[grep("tRNA",annotation.list)]#subset in the parsed annotation
        for (tRNA.term.index in 1:length(unique.tRNA)){
          current.locus<-subset(tRNA.Search, tRNA.Search$Locus==unique.tRNA[tRNA.term.index])#subset the tRNA terms by the current locus
          synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
          if (unique.tRNA[tRNA.term.index] %in% Duplicates ==TRUE){#if the locus is a duplicate
            current.dup<-subset(dup.frame, dup.frame$Duplicates==as.character(unique.tRNA[tRNA.term.index]))
            for (synonym.index in 1:length(synonyms)){
              found.tRNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)
              if (length(found.tRNA)>0){
                max.instance<-min(c(length(found.tRNA),current.dup$DuplicateInstances))#to control number found and written, get lowest common number
                for (dup.found.index in 1:max.instance){
                  found.seq<-seqinr::getSequence(rec$req[[found.tRNA[dup.found.index]]])####Trans work?
                  seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.tRNA[tRNA.term.index],dup.found.index,".fasta"),open="a")
                  Accession.Table[accession.index,grep(paste0("\\b",unique.tRNA[tRNA.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                }
                break
              }
            }
          }
          else{for (synonym.index in 1:length(synonyms)){
            found.tRNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
            if (length(found.tRNA)>0){
              found.seq<-seqinr::getSequence(rec$req[found.tRNA[1]], as.string=FALSE)
              seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.tRNA[tRNA.term.index],".fasta"),open="a")
              Accession.Table[accession.index,grep(paste0("\\b",unique.tRNA[tRNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
              break}}
          }
        }
      }
      if (uni.type[loci.type.index]=="CDS")  {
        rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=CDS", sep=" "))#get the CDS
        current.annot<-annotation.list[grep("CDS",annotation.list)]#subset in the parsed annotation
        for (CDS.term.index in 1:length(unique.CDS)){
          current.locus<-subset(CDS.Search, CDS.Search$Locus==unique.CDS[CDS.term.index])#subset the CDS terms by the current locus
          synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
          if (unique.CDS[CDS.term.index] %in% Duplicates ==TRUE){#if the locus is a duplicate
            current.dup<-subset(dup.frame, dup.frame$Duplicates==as.character(unique.CDS[CDS.term.index]))
            for (synonym.index in 1:length(synonyms)){
              found.CDS<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)
              if (length(found.CDS)>0){
                max.instance<-min(c(length(found.CDS),current.dup$DuplicateInstances))#to control number found and written, get lowest common number
                for (dup.found.index in 1:max.instance){
                  if (TranslateSeqs==TRUE){
                    found.seq<-seqinr::getTrans(rec$req[[found.CDS[dup.found.index]]], numcode=TranslateCode)####Trans work?
                  }
                  else{found.seq<-seqinr::getSequence(rec$req[found.CDS[1]], as.string=FALSE)}
                  seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.CDS[CDS.term.index],dup.found.index,".fasta"),open="a")
                  Accession.Table[accession.index,grep(paste0("\\b",unique.CDS[CDS.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                }
                break
              }
            }
          }
          else{for (synonym.index in 1:length(synonyms)){
            found.CDS<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
            if (length(found.CDS)>0){
              ifelse(TranslateSeqs==TRUE, found.seq<-seqinr::getTrans(rec$req[[found.CDS]],numcode=TranslateCode),found.seq<-seqinr::getSequence(rec$req[found.CDS[1]], as.string=FALSE))
              seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.CDS[CDS.term.index],".fasta"),open="a")
              Accession.Table[accession.index,grep(paste0("\\b",unique.CDS[CDS.term.index],"\\b"), colnames(Accession.Table))]<-new.access
              break}}
          }
        }
      }
      if (uni.type[loci.type.index]=="rRNA")  {
        rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=rRNA", sep=" "))#get the rRNA
        current.annot<-annotation.list[grep("rRNA",annotation.list)]#subset in the parsed annotation
        for (rRNA.term.index in 1:length(unique.rRNA)){
          current.locus<-subset(rRNA.Search, rRNA.Search$Locus==unique.rRNA[rRNA.term.index])#subset the rRNA terms by the current locus
          synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
          if (unique.rRNA[rRNA.term.index] %in% Duplicates ==TRUE){#if the locus is a duplicate
            current.dup<-subset(dup.frame, dup.frame$Duplicates==as.character(unique.rRNA[rRNA.term.index]))
            for (synonym.index in 1:length(synonyms)){
              found.rRNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)
              if (length(found.rRNA)>0){
                max.instance<-min(c(length(found.rRNA),current.dup$DuplicateInstances))#to control number found and written, get lowest common number
                for (dup.found.index in 1:max.instance){
                  found.seq<-seqinr::getSequence(rec$req[[found.rRNA[dup.found.index]]])####Trans work?
                  seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.rRNA[rRNA.term.index],dup.found.index,".fasta"),open="a")
                  Accession.Table[accession.index,grep(paste0("\\b",unique.rRNA[rRNA.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                }
                break
              }
            }
          }
          else{for (synonym.index in 1:length(synonyms)){
            found.rRNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
            if (length(found.rRNA)>0){
              found.seq<-seqinr::getSequence(rec$req[found.rRNA[1]], as.string=FALSE)
              seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.rRNA[rRNA.term.index],".fasta"),open="a")
              Accession.Table[accession.index,grep(paste0("\\b",unique.rRNA[rRNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
              break}}
          }
        }
      }
      if (uni.type[loci.type.index]=="misc_RNA")  {
        rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=misc_RNA", sep=" "))#get the misc_RNA
        current.annot<-annotation.list[grep("misc_RNA",annotation.list)]#subset in the parsed annotation
        for (misc_RNA.term.index in 1:length(unique.misc_RNA)){
          current.locus<-subset(misc_RNA.Search, misc_RNA.Search$Locus==unique.misc_RNA[misc_RNA.term.index])#subset the misc_RNA terms by the current locus
          synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
          if (unique.misc_RNA[misc_RNA.term.index] %in% Duplicates ==TRUE){#if the locus is a duplicate
            current.dup<-subset(dup.frame, dup.frame$Duplicates==as.character(unique.misc_RNA[misc_RNA.term.index]))
            for (synonym.index in 1:length(synonyms)){
              found.misc_RNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)
              if (length(found.misc_RNA)>0){
                max.instance<-min(c(length(found.misc_RNA),current.dup$DuplicateInstances))#to control number found and written, get lowest common number
                for (dup.found.index in 1:max.instance){
                  found.seq<-seqinr::getSequence(rec$req[[found.misc_RNA[dup.found.index]]])####Trans work?
                  seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.misc_RNA[misc_RNA.term.index],dup.found.index,".fasta"),open="a")
                  Accession.Table[accession.index,grep(paste0("\\b",unique.misc_RNA[misc_RNA.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                }
                break
              }
            }
          }
          else{for (synonym.index in 1:length(synonyms)){
            found.misc_RNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
            if (length(found.misc_RNA)>0){
              found.seq<-seqinr::getSequence(rec$req[found.misc_RNA[1]], as.string=FALSE)
              seqinr::write.fasta(found.seq,names=seq.name, paste0(unique.misc_RNA[misc_RNA.term.index],".fasta"),open="a")
              Accession.Table[accession.index,grep(paste0("\\b",unique.misc_RNA[misc_RNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
              break}}
          }
        }
      }
      if (uni.type[loci.type.index]=="D-loop")  {
        mito.loop <- seqinr::query("mito.loop",paste0("AC=",new.access), virtual = TRUE)
        dloop <- seqinr::extractseqs("mito.loop", operation = "feature", feature = "D-loop")
        if (length(dloop)>0){
          dloop.fasta <- seqinr::read.fasta(textConnection(dloop))
          seqinr::write.fasta(dloop.fasta,file="D_loop.fasta",names=seq.name, open="a")
          Accession.Table[accession.index,grep(paste0("\\b","D_loop","\\b"), colnames(Accession.Table))]<-new.access
        }
      }
    }
  }
  #Make Final Accession Table
  UniqueSpecies<-unique(Accession.Table$Species)
  Final.Accession.Table<-data.frame(data.frame(matrix(NA, nrow = length(UniqueSpecies), ncol = 1+length(file.names))))#make empty accession table
  colnames(Final.Accession.Table)<-c("Species",file.names)#attach loci names as colnames
  for (species.index in 1:length(UniqueSpecies)){
    current.spec<-subset(Accession.Table,Accession.Table$Species==UniqueSpecies[species.index])
    Final.Accession.Table[species.index,1]<-UniqueSpecies[species.index]
    for (gene.index in 2:length(Accession.Table)){
      current.loci<-current.spec[,gene.index]#The current column/loci
      found.accessions<-subset(current.loci,!is.na(current.loci))#subset the ones that are present
      numbers<-ifelse(length(found.accessions)==0, NA, paste(found.accessions, sep=",", collapse=","))#ifelse statemen, if not present, fill with NA
      Final.Accession.Table[species.index,gene.index]<-numbers
    }
  }
  Sort.Final.Accession.Table<-Final.Accession.Table[order(Final.Accession.Table$Species),]
  Sort.Final.Accession.Table
}
