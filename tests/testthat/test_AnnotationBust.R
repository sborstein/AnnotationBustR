#Test that the main function AnnotationBust is in fact working properly
#Run on one mitogenome which is missing ND5, which should return NA
#Run with duplicate genes and translated seqs
test_that("AnnotationBust Works",{
skip_on_cran()
data(mtDNAterms)
#load("sysdata.rda")#Manual output of all genes. Missing ND5 gene (misc_feature). Should be missing in real run.
#test with accession=JN628859.1

AnnotationBustOut <- structure(list(Species = "Pallidochromis_tokolosh", tRNA_Phe = "JN628859.1",
    rRNA_12S = "JN628859.1", tRNA_Val = "JN628859.1", rRNA_16S = "JN628859.1",
    ND1 = "JN628859.1", tRNA_Ile = "JN628859.1", tRNA_Gln = "JN628859.1",
    tRNA_Met = "JN628859.1", ND2 = "JN628859.1", tRNA_Trp = "JN628859.1",
    tRNA_Ala = "JN628859.1", tRNA_Asn = "JN628859.1", tRNA_Cys = "JN628859.1",
    tRNA_Tyr = "JN628859.1", COI = "JN628859.1", tRNA_Asp = "JN628859.1",
    COII = "JN628859.1", tRNA_Lys = "JN628859.1", ATP8 = "JN628859.1",
    ATP6 = "JN628859.1", COIII = "JN628859.1", tRNA_Gly = "JN628859.1",
    ND3 = "JN628859.1", tRNA_Arg = "JN628859.1", ND4L = "JN628859.1",
    ND4 = "JN628859.1", tRNA_His = "JN628859.1", ND5 = NA, ND6 = "JN628859.1",
    tRNA_Glu = "JN628859.1", CYTB = "JN628859.1", tRNA_Thr = "JN628859.1",
    tRNA_Pro = "JN628859.1", D_loop = "JN628859.1", tRNA_Ser1 = "JN628859.1",
    tRNA_Ser2 = "JN628859.1", tRNA_Leu1 = "JN628859.1", tRNA_Leu2 = "JN628859.1",Reference = "Darrin Hulsey,C., Keck,B.P., Alamillo,H. and O'Meara,B.C. Mitochondrial genome primers for Lake Malawi cichlids. Mol Ecol Resour 13 (3), 347-353 (2013)."), names = c("Species",
"tRNA_Phe", "rRNA_12S", "tRNA_Val", "rRNA_16S", "ND1", "tRNA_Ile",
"tRNA_Gln", "tRNA_Met", "ND2", "tRNA_Trp", "tRNA_Ala", "tRNA_Asn",
"tRNA_Cys", "tRNA_Tyr", "COI", "tRNA_Asp", "COII", "tRNA_Lys",
"ATP8", "ATP6", "COIII", "tRNA_Gly", "ND3", "tRNA_Arg", "ND4L",
"ND4", "tRNA_His", "ND5", "ND6", "tRNA_Glu", "CYTB", "tRNA_Thr",
"tRNA_Pro", "D_loop", "tRNA_Ser_1", "tRNA_Ser_2", "tRNA_Leu_1",
"tRNA_Leu_2", "Reference"), class = "data.frame", row.names = c(NA, -1L))
setwd (tempdir())
test.out<-try(AnnotationBust(Accessions = "JN628859.1", Terms = mtDNAterms, Duplicates = c("tRNA_Leu","tRNA_Ser"), DuplicateInstances = c(2,2), TranslateSeqs = "Both"))
if ("try-error"%in%class(test.out)) {
    skip("could not connect to remote database")
}else{
    expect_identical(AnnotationBustOut[1,],test.out[1,])
}
#Test it dies if duplicates are off.
expect_error(AnnotationBust(Accessions = "JN628859.1", Terms = mtDNAterms, Duplicates = c("tRNA_Leu","tRNA_Ser"), DuplicateInstances = c(2,2,2), TranslateSeqs = "Only"), "Length of Duplicates and DuplicateInstances is not equal. You must specify the number of duplicates for each feature with duplicates you would like to extract")

#Test no duplicates Features and untidy accession table, and prefix
ncbi.accessions <- c("FJ706295","FJ706343","FJ706292")
data(rDNAterms)#load rDNA search terms from AnnotationBustR
my.sequences <- AnnotationBust(Accessions = ncbi.accessions, rDNAterms, DuplicateSpecies=TRUE, 
Prefix="Example1", Reference = TRUE, TidyAccessions = FALSE)
expect_equal(nrow(my.sequences), 3)#Check number of rows
expect_equal("FJ706343", my.sequences$ITS1[2])#Check number of rows

#Test internet connectivity as well as bad accessions
expect_message(httptest::without_internet(AnnotationBust(Accessions = ncbi.accessions, rDNAterms, DuplicateSpecies=TRUE, 
                                                       Prefix="Example1", Reference = TRUE, TidyAccessions = FALSE)), regexp = "Could not connect to NCBI")

expect_message(AnnotationBust(Accessions = "ZZ999999", rDNAterms, DuplicateSpecies=TRUE, 
                                                         Prefix="Example1", Reference = TRUE, TidyAccessions = FALSE))

#####Check for Introns/Exons####
#Subset out matK from cpDNAterms
cds.terms <- subset(cpDNAterms,cpDNAterms$Feature=="matK")
#Create a vecotr of NA so we can merge with the search terms for introns and exons
cds.terms <- cbind(cds.terms,(rep(NA,length(cds.terms$Feature))))
colnames(cds.terms)[4] <- "IntronExonNumber"

IntronExon.terms<-subset(cpDNAterms,cpDNAterms$Feature=="trnK")
 
##As we want to go for two exons, we will want the synonyms repeated as we are doing and intron and an exon
IntronExon.terms<-rbind(IntronExon.terms,IntronExon.terms)#duplicate the terms
#rep the sequence type we want to extract
IntronExon.terms$Type <- rep(c("intron","intron","exon","exon"))
IntronExon.terms$Feature <- rep(c("trnK_Intron","trnK_Exon2"),each=2)
IntronExon.terms <- cbind(IntronExon.terms,rep(c(1,1,2,2)))#Add intron/exon number info
 
#change column name for number info for IntronExon name
colnames(IntronExon.terms)[4] <- "IntronExonNumber"
 
#We can then merge everything together with MergeSearchTerms terms
IntronExonExampleTerms <- MergeSearchTerms(IntronExon.terms,cds.terms)

IntronExon.example <- AnnotationBust(Accessions=c("KX687910.1"),
Terms=IntronExonExampleTerms, Prefix="DemoIntronExon")

expect_all_equal(unlist(IntronExon.example[,2:4]),"KX687910.1")

####CHECK FEWER DUPLICATES THAN SPECIFIED####
expect_warning(AnnotationBust(Accessions = "JN628859.1", Terms = mtDNAterms, Duplicates = c("tRNA_Leu","tRNA_Ser"), DuplicateInstances = c(2,3), TranslateSeqs = "Both"))

####Test Duplicated CDS and translations####
DupCDS_ex <- AnnotationBust(Accessions = "PX244600.1", Terms = cpDNAterms[cpDNAterms$Name=="ndhB",], TranslateSeqs = "Both", Duplicates = c("ndhB"), DuplicateInstances = c(2))
expect_all_equal(unlist(DupCDS_ex[,2:3]),"PX244600.1")


}

)
