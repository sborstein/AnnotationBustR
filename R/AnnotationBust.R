#' Breaks up genbank sequences into their annotated components based on a set of search terms and writes each subsequence of interest to a FASTA for each accession number supplied.
#' @param Accessions A vector of INSDC GenBank accession numbers. Note that refseq numbers are not supported at the moment (i.e. prefixes like XM_ and NC_) will not work.
#' @param Terms A data frame of search terms. Search terms for animal mitogenomes, nuclear rRNA, chloroplast genomes, and plant mitogenomes are pre-made and can be loaded using the data()function. Additional terms can be added using the MergeSearchTerms function.
#' @param Duplicates A vector of loci names that have more than one copy. Default is NULL
#' @param DuplicateInstances A numeric vector the length of Duplicates of the number of duplicates for each duplicated subsequence you wish to extract.Default is NULL
#' @param TranslateSeqs Should the sequence should be translated to the corresponding peptide sequence. Options include, TRUE, FALSE, or Both. Both returns both the DNA sequence and the corresponding amino acids. Default is FALSE. Note that this only makes sense to list as TRUE or Both for protein coding sequences.
#' @param DuplicateSpecies Logical. As to whether there are duplicate individuals per species. If TRUE, adds the accession number to the fasta file
#' @param Prefix Character If a prefix is specified, all output FASTA files written will begin with the prefix. Default is NULL.
#' @param TidyAccessions Logical as to whether the accession table should have a single row per species. If numerous accessions for a species occur, they will be separated by a comma in the accession table. Default=TRUE.
#' @details The AnnotationBust function takes a vector of accession numbers and a data frame of search terms and extracts subsequences from genomes or concatenated sequences.
#' This function requires a steady internet connection. It writes files in the FASTA format to the working directory and returns an accession table. Files append, so use different prefixes between runs, otherwise they will be added to the current files in the working directory if the names are the same. AnnoitationBustR comes with pre-made
#' search terms for metazoan mitogenomes, plant mitogenomes, chloroplast genomes, and rDNA that can be loaded using data(mtDNAterms), data(mtDNAtermsPlants),data(cpDNAterms), and data(rDNAterms) respectively.
#' Search terms can be completely made by the user as long as they follow a similar format with three columns. The first, Locus, should contain the name of the locus that will also be used to name the files. We recommend following
#' a similar naming convention to what we currently have in the pre-made data frames to ensure that files are named properly, characters like "-" or ".", and names starting with numbers should be avoided as it can cause errors with R.
#' The second column, Type, contains the type of subsequence it is, with options being CDS, rRNA, tRNA, misc_RNA, Intro, Exon, and D_Loop. The last column, Name, consists of a possible
#' name for the locus of interest as it might appear in an annotation. For numerous synonyms for the same locus, one should have each synonym as its own row.An additional fourth column is needed for extracting introns/exons.
#' This column, called IntronExonNumber should contain the number of the intron to extract. If extracting both introns/exons and non-intron/exon sequences the fourth column should be NA for non-intron/exon sequence types. See the examples below and the vignette for detailed examples on extracting intron and exons.
#' It is possible that some subsequences are not fully annotated on ACNUC and, therefore, are not extractable. These will return in the accession table as "type not fully Ann". It is also possible that the sequence has no annotations at all, for which it will return "No Ann. For".
#' If the function returns "Acc. Not Found", the accession number supplied could not be found on NCBI. If "Not On ACNUC GenBank" is returned, the accession is not available through AcNUC.
#' This may be due to ACNUC not being fully up to date. To see the last time ACNUC was updated, run seqinr::choosebank("genbank", infobank=T).
#'
#' For a more detailed walkthrough on using AnnotationBust you can access the vignette with vignette("AnnotationBustR).
#' @return Writes a fasta file(s) to the current working directory selected for each unique subsequence of interest in Terms containing all the accession numbers the subsequence was found in.
#' @return An accesion table of class data.frame.
#' @references Borstein, Samuel R., and Brian C. O'Meara. "AnnotationBustR: An R package to extract subsequences from GenBank annotations." PeerJ Preprints 5 (2017): e2920v1.
#' @examples
#' \dontrun{
#' #Create vector of three NCBI accessions of rDNA toget subsequences of and load rDNA terms.
#' ncbi.accessions<-c("FJ706295","FJ706343","FJ706292")
#' data(rDNAterms)#load rDNA search terms from AnnotationBustR
#' #Run AnnotationBustR and write files to working directory
#' my.sequences<-AnnotationBust(ncbi.accessions, rDNAterms, DuplicateSpecies=TRUE,
#' Prefix="Example1")
#' my.sequences#Return the accession table for each species.
#'
#' ###Example With matK CDS and trnK Introns/Exons##
#' #Subset out matK from cpDNAterms
#' cds.terms<-subset(cpDNAterms,cpDNAterms$Locus=="matK")
#' #Create a vecotr of NA so we can merge with the search terms for introns and exons
#' cds.terms<-cbind(cds.terms,(rep(NA,length(cds.terms$Locus))))
#' colnames(cds.terms)[4]<-"IntronExonNumber"
#'
#' #Prepare a search term table for the intron and exons to remove
#' #We can start with the cpDNAterms for trnK
#' IntronExon.terms<-subset(cpDNAterms,cpDNAterms$Locus=="trnK")
#' #As we want to go for two exons, we will want the synonyms repeated as we are doing and intron
#' #and an exon
#' IntronExon.terms<-rbind(IntronExon.terms,IntronExon.terms)#duplicate the terms
#' #rep the sequence type we want to extract
#' IntronExon.terms$Type<-rep(c("Intron","Intron","Exon","Exon"))
#' IntronExon.terms$Locus<-rep(c("trnK_Intron","trnK_Exon2"),each=2)
#' IntronExon.terms<-cbind(IntronExon.terms,rep(c(1,1,2,2)))#Add intron/exon number info
#' #change column name for number info for IntronExon name
#' colnames(IntronExon.terms)[4]<-"IntronExonNumber"
#' #We can then merge everything together with MergeSearchTerms terms
#' IntronExonExampleTerms<-MergeSearchTerms(IntronExon.terms,cds.terms)
#'
#' #Run AnnotationBust
#' IntronExon.example<-AnnotationBust(Accessions=c("KX687911.1", "KX687910.1"),
#' Terms=IntronExonExampleTerms, Prefix="DemoIntronExon")
#' }
#' @export

AnnotationBust<-function(Accessions, Terms, Duplicates= NULL,DuplicateInstances=NULL, TranslateSeqs=FALSE, TranslateCode=1, DuplicateSpecies=FALSE, Prefix=NULL, TidyAccessions=TRUE){
  #seqinr::choosebank("genbank")
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
  #Setup Prefix
  if(is.null(Prefix)){
    File.Prefix<-NULL
  }else{File.Prefix<-paste0(Prefix,"_")}

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
      misc_RNA.Search<-subset(Terms, Terms$Type=="misc_RNA")
      unique.misc_RNA<-unique(misc_RNA.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="Intron"){
      Intron.Search<-subset(Terms, Terms$Type=="Intron")
      unique.Intron<-unique(Intron.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="Exon"){
      Exon.Search<-subset(Terms, Terms$Type=="Exon")
      unique.Exon<-unique(Exon.Search$Locus)
    }
    if (uni.type[subsequence.type.index]=="misc_feature"){
      misc_feature.Search<-subset(Terms, Terms$Type=="misc_feature")
      unique.misc_feature<-unique(misc_feature.Search$Locus)
    }
  }
  for (accession.index in 1:length(Accessions)){
    try.test<-try(seqinr::choosebank("genbank", verbose = FALSE),silent = TRUE)
    while(class(try.test)== "try-error"){
      Sys.sleep(10)
      try.test<-try(seqinr::choosebank("genbank", verbose = FALSE),silent = TRUE)
    }
    new.access<-strsplit(Accessions[accession.index],"\\.",perl=TRUE)[[1]][1]#split and decimal spot in accession number. seqinr won't take them with it
    species.name<-try(attr(ape::read.GenBank(Accessions[accession.index]),"species"))#get the sequence names
    if(class(species.name)=="try-error"){
      Accession.Table[accession.index,1:length(colnames(Accession.Table))]<-paste(new.access,"Acc. Not Found",sep=" ")
      next
    }else{
      print(paste("Working On Accession",new.access,species.name, sep=" "))
      ifelse(DuplicateSpecies==TRUE, seq.name<-paste(species.name,new.access,sep = "_"), seq.name<-species.name)
      Accession.Table$Species[accession.index]<-species.name
      full.rec<-try(seqinr::query(paste0("AC=",new.access)))
      if(grepl("unknown accession", full.rec[[1]][1])==TRUE){
        Accession.Table[accession.index,2:length(colnames(Accession.Table))]<-paste(new.access,"Not On ACNUC GenBank",sep=" ")
        next
      }else{
        while(class(full.rec)== "try-error"){
          Sys.sleep(10)
          try.test<-try(seqinr::choosebank("genbank", verbose = FALSE),silent = TRUE)
          full.rec<-try(seqinr::query(paste0("AC=",new.access)))
        }
      }
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
      #check that annotation is multiple parts
      check.ann<-annotation.list[-1]
      if(length(check.ann)==1){
        Accession.Table[accession.index,2:length(colnames(Accession.Table))]<-paste("No Ann. For",new.access,sep=" ")
        next
      }
      #Now find loci for every loci type
      for (loci.type.index in 1:length(uni.type)){
        if (uni.type[loci.type.index]=="tRNA")  {
          rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=tRNA", sep=" "))#get the tRNA
          current.annot<-annotation.list[grep("tRNA            ",annotation.list)]#subset in the parsed annotation
          if (!length(current.annot)==length(rec$req)){
            bad.cols<-which(colnames(Accession.Table) %in% unique.tRNA)
            Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
            next
          }
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
                    found.seq<-seqinr::getSequence(rec$req[[found.tRNA[dup.found.index]]])#
                    seqinr::write.fasta(toupper(found.seq),names=seq.name, paste0(File.Prefix,unique.tRNA[tRNA.term.index],dup.found.index,".fasta"),open="a", nbchar = 60)
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
                seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.tRNA[tRNA.term.index],".fasta"),open="a",nbchar = 60)
                Accession.Table[accession.index,grep(paste0("\\b",unique.tRNA[tRNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                break}}
            }
          }
        }
        if (uni.type[loci.type.index]=="CDS")  {
          rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=CDS", sep=" "))#get the CDS
          current.annot<-annotation.list[grep("CDS            ",annotation.list)]#subset in the parsed annotation
          if (!length(current.annot)==length(rec$req)){
            bad.cols<-which(colnames(Accession.Table) %in% unique.CDS)
            Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
            next
          }
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
                      found.seq<-seqinr::c2s(seqinr::getTrans(rec$req[[found.CDS[dup.found.index]]], numcode=TranslateCode))####Trans work?
                    }
                    else{found.seq<-seqinr::getSequence(rec$req[found.CDS[dup.found.index]], as.string=FALSE)}
                    seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.CDS[CDS.term.index],dup.found.index,".fasta"),open="a",nbchar = 60)
                    Accession.Table[accession.index,grep(paste0("\\b",unique.CDS[CDS.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                  }
                  break
                }
              }
            }
            else{for (synonym.index in 1:length(synonyms)){
              found.CDS<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
              if (length(found.CDS)>0){
                ifelse(TranslateSeqs==TRUE, found.seq<-seqinr::c2s(seqinr::getTrans(rec$req[[found.CDS]],numcode=TranslateCode)),found.seq<-seqinr::getSequence(rec$req[found.CDS[1]], as.string=FALSE))
                seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.CDS[CDS.term.index],".fasta"),open="a",nbchar = 60)
                Accession.Table[accession.index,grep(paste0("\\b",unique.CDS[CDS.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                break}}
            }
          }
        }
        if (uni.type[loci.type.index]=="rRNA")  {
          rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=rRNA", sep=" "))#get the rRNA
          current.annot<-annotation.list[grep("rRNA            ",annotation.list)]#subset in the parsed annotation
          if (!length(current.annot)==length(rec$req)){
            bad.cols<-which(colnames(Accession.Table) %in% unique.rRNA)
            Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
            next
          }
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
                    seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.rRNA[rRNA.term.index],dup.found.index,".fasta"),open="a",nbchar = 60)
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
                seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.rRNA[rRNA.term.index],".fasta"),open="a",nbchar = 60)
                Accession.Table[accession.index,grep(paste0("\\b",unique.rRNA[rRNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                break}}
            }
          }
        }
        if (uni.type[loci.type.index]=="misc_RNA")  {
          rec<-seqinr::query(paste("SUB", paste0("AC=",new.access), "AND T=misc_RNA", sep=" "))#get the misc_RNA
          current.annot<-annotation.list[grep("misc_RNA        ",annotation.list)]#subset in the parsed annotation
          if (!length(current.annot)==length(rec$req)){
            bad.cols<-which(colnames(Accession.Table) %in% unique.misc_RNA)
            Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
            next
          }
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
                    seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.misc_RNA[misc_RNA.term.index],dup.found.index,".fasta"),open="a",nbchar = 60)
                    Accession.Table[accession.index,grep(paste0("\\b",unique.misc_RNA[misc_RNA.term.index],dup.found.index,"\\b"), colnames(Accession.Table))]<-new.access
                  }
                  break
                }
              }
            }
            else{for (synonym.index in 1:length(synonyms)){
              found.misc_RNA<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#search for the regular
              if (length(found.misc_RNA)>0){
                found.seq<-(seqinr::getSequence(rec$req[found.misc_RNA[1]], as.string=FALSE))
                seqinr::write.fasta(toupper(found.seq[[1]]),names=seq.name, paste0(File.Prefix,unique.misc_RNA[misc_RNA.term.index],".fasta"),open="a", nbchar = 60)
                Accession.Table[accession.index,grep(paste0("\\b",unique.misc_RNA[misc_RNA.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                break}}
            }
          }
        }
        if (uni.type[loci.type.index]=="Intron")  {
          current.annot<-annotation.list[grep("intron          ",annotation.list)]#subset in the parsed annotation
          intron.find <- seqinr::query("intron.find",paste0("AC=",new.access), virtual = TRUE)
          check.intron <- seqinr::extractseqs("intron.find", operation = "feature", feature = "intron")
          if (length(check.intron)>0){
            intron.fasta <- seqinr::read.fasta(textConnection(check.intron),as.string = FALSE,forceDNAtolower = FALSE)
            if (!length(current.annot)==length(intron.fasta)){
              bad.cols<-which(colnames(Accession.Table) %in% unique.Intron)
              Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
              next
            }
            for (Intron.term.index in 1:length(unique.Intron)){
              current.locus<-subset(Intron.Search, Intron.Search$Locus==unique.Intron[Intron.term.index])#subset the Intron terms by the current locus
              synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
              IntronNumber<-unique(current.locus$IntronExonNumber)
              for (synonym.index in 1:length(synonyms)){
                found.Intron<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#Find introns for gene name
                foundNumber<-grep(paste0("\\b", IntronNumber[synonym.index],"\\b"),current.annot[found.Intron])#find intron number that is to be extracted
                ifelse(length(foundNumber)==0|length(found.Intron)==1,foundNumber<-found.Intron,foundNumber<-foundNumber)
                if (length(foundNumber)>0){
                  seqinr::write.fasta(intron.fasta[foundNumber],names=seq.name, paste0(File.Prefix,unique.Intron[Intron.term.index],".fasta"),open="a",nbchar = 60)
                  Accession.Table[accession.index,grep(paste0("\\b",unique.Intron[Intron.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                  break}
              }
            }
          }
        }
        if (uni.type[loci.type.index]=="Exon")  {
          current.annot<-annotation.list[grep("exon            ",annotation.list)]#subset in the parsed annotation
          Exon.find <- seqinr::query("Exon.find",paste0("AC=",new.access), virtual = TRUE)
          check.Exon <- seqinr::extractseqs("Exon.find", operation = "feature", feature = "Exon")
          if (length(check.Exon)>0){
            Exon.fasta <- seqinr::read.fasta(textConnection(check.Exon),as.string = FALSE,forceDNAtolower = FALSE)
            if (!length(current.annot)==length(Exon.fasta)){
              bad.cols<-which(colnames(Accession.Table) %in% unique.Exon)
              Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
              next
            }
            for (Exon.term.index in 1:length(unique.Exon)){
              current.locus<-subset(Exon.Search, Exon.Search$Locus==unique.Exon[Exon.term.index])#subset the Exon terms by the current locus
              synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
              ExonNumber<-unique(current.locus$IntronExonNumber)
              for (synonym.index in 1:length(synonyms)){
                found.Exon<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#Find Exons for gene name
                foundNumber<-grep(paste0("\\b", ExonNumber[synonym.index],"\\b"),current.annot[found.Exon])#find Exon number that is to be extracted
                ifelse(length(foundNumber)==0|length(found.Exon)==1,foundNumber<-found.Exon,foundNumber<-foundNumber)
                if (length(foundNumber)>0){
                  seqinr::write.fasta(Exon.fasta[foundNumber],names=seq.name, paste0(File.Prefix,unique.Exon[Exon.term.index],".fasta"),open="a",nbchar = 60)
                  Accession.Table[accession.index,grep(paste0("\\b",unique.Exon[Exon.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                  break}
              }
            }
          }
        }
        if (uni.type[loci.type.index]=="misc_feature")  {
          current.annot<-annotation.list[grep("misc_feature    ",annotation.list)]#subset in the parsed annotation
          misc_feature.find <- seqinr::query("misc_feature.find",paste0("AC=",new.access), virtual = TRUE)
          check.misc_feature <- seqinr::extractseqs("misc_feature.find", operation = "feature", feature = "misc_feature")
          if (length(check.misc_feature)>0){
            misc_feature.fasta <- seqinr::read.fasta(textConnection(check.misc_feature),as.string = FALSE,forceDNAtolower = FALSE)
            if (!length(current.annot)==length(misc_feature.fasta)){
              bad.cols<-which(colnames(Accession.Table) %in% unique.misc_feature)
              Accession.Table[accession.index,bad.cols]<-paste("type not fully Ann.",new.access,sep=" ")
              next
            }
            for (misc_feature.term.index in 1:length(unique.misc_feature)){
              current.locus<-subset(misc_feature.Search, misc_feature.Search$Locus==unique.misc_feature[misc_feature.term.index])#subset the misc_feature terms by the current locus
              synonyms<-unique(current.locus$Name)#subset the Name column, which includes the synonyms
              for (synonym.index in 1:length(synonyms)){
                found.misc_feature<-grep(paste0("\\b",synonyms[synonym.index],"\\b"), current.annot)#Find Exons for gene name
                if (length(found.misc_feature)>0){
                  seqinr::write.fasta(misc_feature.fasta[found.misc_feature],names=seq.name, paste0(File.Prefix,unique.misc_feature[misc_feature.term.index],".fasta"),open="a",nbchar = 60)
                  Accession.Table[accession.index,grep(paste0("\\b",unique.misc_feature[misc_feature.term.index],"\\b"), colnames(Accession.Table))]<-new.access
                  break}
              }
            }
          }
        }
        if (uni.type[loci.type.index]=="D-loop")  {
          mito.loop <- seqinr::query("mito.loop",paste0("AC=",new.access), virtual = TRUE)
          dloop <- seqinr::extractseqs("mito.loop", operation = "feature", feature = "D-loop")
          if (length(dloop)>0){
            dloop.fasta <- seqinr::read.fasta(textConnection(dloop),as.string = FALSE,forceDNAtolower = FALSE)
            seqinr::write.fasta(dloop.fasta,file=paste0(File.Prefix,"D_loop.fasta"),names=seq.name, open="a")
            Accession.Table[accession.index,grep(paste0("\\b","D_loop","\\b"), colnames(Accession.Table))]<-new.access
          }
        }
      }
      seqinr::closebank()
    }
  }
  #Make Final Accession Table
  if(TidyAccessions==TRUE){
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
      Sort.Final.Accession.Table<-Final.Accession.Table[order(Final.Accession.Table$Species),]
    }
    rownames(Sort.Final.Accession.Table)<-1:nrow(Sort.Final.Accession.Table)
  }
  }else{
  Final.Accession.Table<-Accession.Table
  Sort.Final.Accession.Table<-Final.Accession.Table[order(Final.Accession.Table$Species),]
  rownames(Sort.Final.Accession.Table)<-1:nrow(Sort.Final.Accession.Table)

}
  Sort.Final.Accession.Table
}




seqinr_getSequence_loop <- function(x) {
  seq_return <- NA
  reps <- 0
  while(is.na(seq_return[1])) {
    reps <- reps + 1
    try(seq_return <- seqinr::getSequence(x))
    if(is.na(seq_return[1])) {
      warning(paste0("Connection to the database failed; trying again in ", 10*reps, " seconds"), immediate.=TRUE)
      Sys.sleep(10*reps)

    }
  }
  return(seq_return)
}
