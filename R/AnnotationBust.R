#' Breaks up genbank sequences into their annotated components based on a set of search terms and writes each subsequence of interest to a FASTA for each accession number supplied.
#' @param Accessions A vector of GenBank accession numbers.
#' @param Terms A data frame of search terms. Search terms for animal mitogenomes, nuclear rRNA, chloroplast genomes, and plant mitogenomes are pre-made and can be loaded using the data() function. Additional terms can be added using the MergeSearchTerms function or a user supplied data frame may be provided.
#' @param Duplicates A vector of the features  which occur more than once in the sequence and for which the extraction of multiple copies is desired. Default is NULL. In the case more than one copy exists when set to NULL, only the first instance will be extracted.
#' @param DuplicateInstances A numeric vector the length of Duplicates of the number of duplicates for each duplicated feature you wish to extract. Only needs to be used if Duplicates are provided.Default is NULL (i.e. no duplicates).
#' @param TranslateSeqs Should coding sequences (cds) be translated to the corresponding peptide sequence? Options include, Only, None, or Both. Both returns both the DNA sequence and the corresponding peptide sequence. Default is FALSE.
#' @param DuplicateSpecies Logical. As to whether there are duplicate individuals per species. If TRUE, adds the accession number to the fasta header when writing sequences to file.
#' @param Prefix Character. If a prefix is specified, all output FASTA files written will begin with the prefix followed by an underscore. Default is NULL (i.e. no prefix).
#' @param TidyAccessions Logical. Should the accession table  have a single row per species? If numerous accessions for a species occur, they will be separated by a comma in the accession table. Default=TRUE.
#' @param Verbose Logical. Should progress be printed to the screen. The current accession and species name will be printed to the screen.
#' @details 
#' The AnnotationBust function takes a vector of accession numbers and a data frame of search terms and extracts subsequences from genomes or concatenated sequences. This function connects directly to the NCBI database and requires a steady internet connection. The function writes files in the FASTA format to the current working directory and returns an accession table. Files append, so use different prefixes between runs, otherwise they will be added to current files in the working directory with the same name. 
#' 
#' AnnotationBustR comes with pre-made search terms for metazoan mitogenomes, plant mitogenomes, chloroplast genomes, and rDNA that can be loaded using data(mtDNAterms), data(mtDNAtermsPlants), data(cpDNAterms), and data(rDNAterms) respectively. Search terms can be completely made by the user if they follow a similar format with three columns. The first, Feature, should contain the name of the features as the user wishes it be displayed in their files as it is used to name the files and create the accession table. We recommend following a similar naming convention to what we currently have in the pre-made data frames to ensure that files are named properly, characters like "-" or ".", and names starting with numbers should be avoided as it can cause errors with R. The second column, Type, contains the type of subsequence it is (eg. CDS, exon, intron, rRNA, tRNA, misc_RNA, misc_feature, D_Loop). The last column, Name, consists of a possible synonym for the feature of interest as it might appear in an annotation. For numerous synonyms for the same locus, one should have each synonym as its own row. An additional fourth column is needed for extracting introns/exons. This column, called IntronExonNumber should contain the number of the desired intron or exon to extract. If extracting both introns/exons and non-intron/exon sequences the fourth column should be NA for non-intron/exon sequence types. See the examples below and the vignette for detailed examples on extracting intron and exons. 
#' 
#' In the event that an accession is not found on NCBI, the message “Accession # not found on NBI”.
#' 
#' For a more detailed walk-through on using AnnotationBust you can access the vignette with vignette("AnnotationBustR).
#' @return Writes a fasta file(s) to the current working directory selected for each unique subsequence of interest in Terms containing all the accession numbers the subsequence was found in.
#' @return An accesion table of class data.frame.
#' @examplesIf interactive()
#' \donttest{
#' ncbi.accessions <- c("FJ706295","FJ706343","FJ706292")
#' data(rDNAterms)#load rDNA search terms from AnnotationBustR
#' my.sequences <- AnnotationBust(Accessions = ncbi.accessions, rDNAterms, DuplicateSpecies=TRUE, 
#' Prefix="Example1")
#' my.sequences
#' 
#' ###Example With matK CDS and addint introns/exons for trnK###
#' #Subset out matK from cpDNAterms
#' cds.terms <- subset(cpDNAterms,cpDNAterms$Feature=="matK")
#' #Create a vecotr of NA so we can merge with the search terms for introns and exons
#' cds.terms <- cbind(cds.terms,(rep(NA,length(cds.terms$Feature))))
#' colnames(cds.terms)[4] <- "IntronExonNumber"
#' 
#' #Prepare a search term table for the intron and exons to remove
#' #We can start with the cpDNAterms for trnK
#' IntronExon.terms<-subset(cpDNAterms,cpDNAterms$Feature=="trnK")
#' 
#' #As we want to go for two exons, we will want the synonyms repeated as we are doing and intron
#' #and an exon
#' IntronExon.terms<-rbind(IntronExon.terms,IntronExon.terms)#duplicate the terms
#' 
#' #rep the sequence type we want to extract
#' IntronExon.terms$Type <- rep(c("intron","intron","exon","exon"))
#' IntronExon.terms$Feature <- rep(c("trnK_Intron","trnK_Exon2"),each=2)
#' IntronExon.terms <- cbind(IntronExon.terms,rep(c(1,1,2,2)))#Add intron/exon number info
#' 
#' #change column name for number info for IntronExon name
#' colnames(IntronExon.terms)[4] <- "IntronExonNumber"
#' 
#' #We can then merge everything together with MergeSearchTerms terms
#' IntronExonExampleTerms <- MergeSearchTerms(IntronExon.terms,cds.terms)
#'
#' #Run AnnotationBust
#' IntronExon.example <- AnnotationBust(Accessions=c("KX687911.1", "KX687910.1"),
#' Terms=IntronExonExampleTerms, Prefix="DemoIntronExon")
#' }
#' @references Borstein SR, and O’Meara BC. 2018. AnnotationBustR: an R package to extract subsequences from GenBank annotations. PeerJ 6:e5179. 10.7717/peerj.5179.
#' @author Samuel R. Borstein, Brian C. O'Meara 
#' @export

AnnotationBust <- function(Accessions, Terms, Duplicates = NULL, DuplicateInstances = NULL, 
                               TranslateSeqs = "None", DuplicateSpecies = FALSE,  
                               Prefix = NULL, TidyAccessions = TRUE, Verbose = TRUE){
  
  uni.feature <- unique(Terms$Feature)#Subset unique Loci

  if(is.null(Duplicates)==FALSE){#Routine for naming if DuplicateInstances==TRUE
    if(!length(Duplicates)==length(DuplicateInstances)){#Check and error if Duplicates does not equal DuplicateInstances
      stop("Length of Duplicates and DuplicateInstances is not equal. You must specify the number of duplicates for each feature with duplicates you would like to extract")
    }
    new.file.names <- character(length = 0)#prep file names
    dup.frame <- data.frame(Duplicates, DuplicateInstances)#prep duplicates
    singles <- uni.feature[uni.feature %in% dup.frame$Duplicates==FALSE]#get matching duplicates
    doubles <- uni.feature[uni.feature %in% dup.frame$Duplicates==TRUE]#get matching duplicates
    for (i in 1:length (doubles)){#for each duplicate
      sub.dups <- subset(dup.frame, dup.frame$Duplicates %in% doubles[i]==TRUE)#subset the duplicate genes
      number.names <- paste0(sub.dups$Duplicates,"_", 1:sub.dups$DuplicateInstances)#add duplicate numbers to names
      new.file.names <- append(new.file.names, number.names)#append the new pasted file names
    }
    file.names <- c(as.vector(singles), new.file.names)#combo the names of single and duplicate loci
  }else {file.names <- as.vector(uni.feature)}#if no duplicates, unique feature names are the file names
  Accession.Table <- data.frame(data.frame(matrix(NA, nrow = length(Accessions), ncol = 1+length(file.names))))#make empty accession table
  colnames(Accession.Table) <- c("Species",file.names)#attach loci names as column names
  if(is.null(Prefix)){#prepare prefix naming for files
    File.Prefix <- NULL#If no prefix, set to NULL
  }else{
    File.Prefix <- paste0(Prefix,"_")#Add underscore if prefix exists for tidyness
    }
  #For each accession
  for(accession.index in 1:length(Accessions)){
    Current.Accession <- Accessions[accession.index]
    raw_gb_text <- try(rentrez::entrez_fetch(db = "nuccore", id = Current.Accession, rettype = "gbwithparts", retmode = "text"))#Use entrez to access accession
    if(inherits(raw_gb_text, "try-error")){#Flag if try error
      message("Could not connect to NCBI")#message error and skip to next accession.
      next
    }
    if(grepl(pattern = "Error: F a i l e d  t o  u n d e r s t a n d",raw_gb_text)){#search for errors
      message(paste0("Working On Accession ", accession.index, " of ", length(Accessions),": ", Accessions[accession.index], " not found on NCBI"))#message if not found at populate the Accession Table
      Accession.Table[accession.index,1:length(colnames(Accession.Table))]<-paste("Failed to find", Accessions[accession.index],sep=" ")#Add not found info
      next#skip to next
    }
    gb_lines <- strsplit(raw_gb_text, "\n|\r\n")[[1]]#string split lines
    org_line_index <- grep("^\\s*ORGANISM", gb_lines)#Find organism line
    org_line <- gb_lines[org_line_index]#Subset line grabbed above for organism
    organism_name <- gsub(" ","_",trimws(sub("^\\s*ORGANISM", "", org_line)))#Extract organism name and replace speces with under scores
    if(Verbose == TRUE){#if verbose, print progress
      message(paste("Working On Accession ", accession.index, " of ", length(Accessions),": ",  Accessions[accession.index],", ",organism_name, sep=""))
    }
    Accession.Table[accession.index,"Species"] <- organism_name
    current_record <- parse_genbank(ReadGB = gb_lines, primary_accession = Current.Accession)#Parse the accession and its annotations.
    ifelse(DuplicateSpecies==TRUE, seq.name<-paste(organism_name,Accessions[accession.index],sep = "_"), seq.name<-organism_name)#For duplicate species, append the accesion to end of species name
    for(term.index in seq_along(uni.feature)){#for each term
      synonyms <- Terms[Terms$Feature==uni.feature[term.index],]#subset synonyms
      for(synonym.index in 1:nrow(synonyms)){#for each synonym
        found.type <- base::grep(pattern = paste0("^",synonyms$Type[synonym.index],"$"),x = names(current_record))#find feature type
        if(length(found.type)==1){#If type is found in annotations for the accession
          term.search <- base::grep(pattern = paste0("^",synonyms$Name[synonym.index],"$"),x = current_record[[found.type]]$gene)#Check genes first
          ifelse(length(term.search)>0, term.search <- term.search, term.search <- base::grep(pattern = paste0("^",synonyms$Name[synonym.index],"$"),x = current_record[[found.type]]$product))#if not in genes, check products
          ifelse(length(term.search)>0, term.search <- term.search, term.search <- base::grep(pattern = paste0("^",synonyms$Name[synonym.index],"$"),x = current_record[[found.type]]$note))#if not in genes, check products
          if(length(term.search)==0 && names(current_record)[[found.type]]=="D-loop"){#Facilitate D-loop discovery
            term.search <- base::grep(pattern = paste0("^",synonyms$Name[synonym.index],"$"),x = current_record[[found.type]]$type)
          }
          if(synonyms[synonym.index,]$Type%in%c("exon","intron")){#If the feature is an intron OR an exon
            temp.term <- current_record[[found.type]][term.search,]#Make a temporary reduced dataset to search
            IntronExonHits <- which(temp.term$number==synonyms$IntronExonNumber[synonym.index])
            term.search <- term.search[IntronExonHits]
          }
          if(length(term.search)>0){#If successful search
            FeatureGrab <- current_record[[found.type]][term.search,]#grab the features
            FeatureIndexCheck <- unique(FeatureGrab$feat_index)#Check to see if there is a single unique ID for annotation or several.
            if(uni.feature[term.index] %in% Duplicates == FALSE){#If current feature does not have duplicate instances requested
              Extracted.seq <- extract_ranges_seq(full_sequence = current_record$sequence, target_gr = FeatureGrab[FeatureGrab$feat_index==FeatureIndexCheck[1]])#Extract Seq based on feature ranges
              names(Extracted.seq) <- seq.name#name sequence for species
              if(synonyms$Type[synonym.index] == "CDS"){#If a CDS, use control flow based on users desire for translation
                if(TranslateSeqs%in%c("Only","Both")){#Translate if desired
                  Current.Trans <- unique(FeatureGrab[FeatureGrab$feat_index==FeatureIndexCheck[1]]$translation)#Get the current translation
                  if(length(Current.Trans) == 1){#ensure only single translation across CDS
                    TransString <- Biostrings::AAStringSet(Current.Trans)#create AA stringset object for translation
                    names(TransString) <- seq.name#name sequence for species
                    Biostrings::writeXStringSet(x = TransString, filepath = paste0(File.Prefix,uni.feature[term.index],"_Translation",".fasta"), format = "fasta",append = T)#Write Translation
                    Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                  }
                }
                if(TranslateSeqs%in%c("Both","None")){#Write the Nucleotides if Both or No Translation
                  Biostrings::writeXStringSet(x = Extracted.seq, filepath = paste0(File.Prefix,uni.feature[term.index],".fasta"), format = "fasta",append = T)#Write DNA sequence
                  Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                }
                break
              }else{#if not a CDS
                Biostrings::writeXStringSet(x = Extracted.seq, filepath = paste0(File.Prefix,uni.feature[term.index],".fasta"), format = "fasta",append = T)#Write DNA sequence
                Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                break
              }
            }else{#if duplicates
              DupID <- which(uni.feature[term.index]==Duplicates)#In case order is off, make sure we get the right one.
              CurrentDupTargets <- DuplicateInstances[DupID]#Subset out the current number of instances desired
              if(CurrentDupTargets>length(unique(FeatureGrab$feat_index))){
                CurrentDupTargets <- length(FeatureGrab)
                warning(paste0("Number of duplicates specified is greater than the number in annotations. Readjusting duplicates for ", uni.feature[term.index]," to ", CurrentDupTargets))
              }
              for(DupIndex in 1:CurrentDupTargets){#for each duplicate
                CurrentDup <- FeatureGrab[FeatureGrab$feat_index == FeatureIndexCheck[DupIndex], ]
                Extracted.seq <- extract_ranges_seq(full_sequence = current_record$sequence, target_gr = CurrentDup[CurrentDup$feat_index==FeatureIndexCheck[DupIndex]])#Extract Seq based on feature ranges
                names(Extracted.seq) <- seq.name#name sequence for species
                if(synonyms$Type[synonym.index] == "cds"){#If a CDS, use control flow based on users desire for translation
                  if(TranslateSeqs%in%c("Only","Both")){#Translate if desired
                    Current.Trans <- unique(CurrentDup$translation)#Get the current translation
                    if(length(Current.Trans) == 1){#ensure only single translation across CDS
                      TransString <- Biostrings::AAStringSet(Current.Trans)#create AA stringset object for translation
                      names(TransString) <- seq.name#name sequence for species
                      Biostrings::writeXStringSet(x = TransString, filepath = paste0(File.Prefix,uni.feature[term.index],"_",DupIndex,"_Translation",".fasta"), format = "fasta",append = T)#Write Translation
                      Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"_",DupIndex,"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                    }
                  }
                  if(TranslateSeqs%in%c("Both","None")){#Write the Nucleotides if Both or No Translation
                    Biostrings::writeXStringSet(x = Extracted.seq, filepath = paste0(File.Prefix,uni.feature[term.index],"_",DupIndex,".fasta"), format = "fasta",append = T)#Write DNA sequence
                    Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"_",DupIndex,"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                  }
                }else{#if not a CDS
                  Biostrings::writeXStringSet(x = Extracted.seq, filepath = paste0(File.Prefix,uni.feature[term.index],"_",DupIndex,".fasta"), format = "fasta",append = T)#Write DNA sequence
                  Accession.Table[accession.index,grep(paste0("\\b",uni.feature[term.index],"_",DupIndex,"\\b"), colnames(Accession.Table))] <- Accessions[accession.index]
                }
              }#for loop duplicate ends
              break#break out of current term
            }#duplicate control end
          }#if length search end
        }#if feature type is found end
      }#end of synonym for loop
    }#for each term loop end
  }#Accession For loop end
  #Make Final Accession Table
  if(TidyAccessions==TRUE){
    UniqueSpecies<-unique(Accession.Table$Species)#Get unique species in table
    Final.Accession.Table<-data.frame(data.frame(matrix(NA, nrow = length(UniqueSpecies), ncol = 1+length(file.names))))#make empty accession table
    colnames(Final.Accession.Table)<-c("Species",file.names)#attach loci names as colnames
    for (species.index in 1:length(UniqueSpecies)){#for each species
      current.spec<-subset(Accession.Table,Accession.Table$Species==UniqueSpecies[species.index])#subset species
      Final.Accession.Table[species.index,1]<-UniqueSpecies[species.index]#Assign species name to final table
      for (gene.index in 2:length(Accession.Table)){#for each feature
        current.loci<-current.spec[,gene.index]#The current column/loci
        found.accessions<-subset(current.loci,!is.na(current.loci))#subset the ones that are present
        numbers<-ifelse(length(found.accessions)==0, NA, paste(found.accessions, sep=",", collapse=","))#ifelse statemen, if not present, fill with NA
        Final.Accession.Table[species.index,gene.index]<-numbers#assign to cells in final table
        Sort.Final.Accession.Table<-Final.Accession.Table[order(Final.Accession.Table$Species),]#Sort by species name
      }
      rownames(Sort.Final.Accession.Table)<-1:nrow(Sort.Final.Accession.Table)#Assign rownames
    }
  }else{
    Final.Accession.Table<-Accession.Table#assign the accession table to return prior to sorting
    Sort.Final.Accession.Table<-Final.Accession.Table[order(Final.Accession.Table$Species),]#sort table by species
    rownames(Sort.Final.Accession.Table)<-1:nrow(Sort.Final.Accession.Table)#assign rownames
  } 
  return(Final.Accession.Table)#return accession table
}#Function End


                      