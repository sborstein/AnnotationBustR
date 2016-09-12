#Test that the main function AnnotationBust is in fact working properly
#Run on one mitogenome which is missing ND5, which should return NA
#Run with duplicate genes and translated seqs
test_that("AnnotationBust Works",{
data(mtDNAterms)
#load("sysdata.rda")#Manual output of all genes. Missing ND5 gene (misc_feature). Should be missing in real run.
#test with accession=JN628859.1

AnnotationBustOut <- structure(list(Species = "Pallidochromis_tokolosh", tRNA_Phe = "JN628859",
    rRNA_12S = "JN628859", tRNA_Val = "JN628859", rRNA_16S = "JN628859",
    ND1 = "JN628859", tRNA_Ile = "JN628859", tRNA_Gln = "JN628859",
    tRNA_Met = "JN628859", ND2 = "JN628859", tRNA_Trp = "JN628859",
    tRNA_Ala = "JN628859", tRNA_Asn = "JN628859", tRNA_Cys = "JN628859",
    tRNA_Tyr = "JN628859", COI = "JN628859", tRNA_Asp = "JN628859",
    COII = "JN628859", tRNA_Lys = "JN628859", ATP8 = "JN628859",
    ATP6 = "JN628859", COIII = "JN628859", tRNA_Gly = "JN628859",
    ND3 = "JN628859", tRNA_Arg = "JN628859", ND4L = "JN628859",
    ND4 = "JN628859", tRNA_His = "JN628859", ND5 = NA, ND6 = "JN628859",
    tRNA_Glu = "JN628859", CYTB = "JN628859", tRNA_Thr = "JN628859",
    tRNA_Pro = "JN628859", D_loop = "JN628859", tRNA_Ser1 = "JN628859",
    tRNA_Ser2 = "JN628859", tRNA_Leu1 = "JN628859", tRNA_Leu2 = "JN628859"), .Names = c("Species",
"tRNA_Phe", "rRNA_12S", "tRNA_Val", "rRNA_16S", "ND1", "tRNA_Ile",
"tRNA_Gln", "tRNA_Met", "ND2", "tRNA_Trp", "tRNA_Ala", "tRNA_Asn",
"tRNA_Cys", "tRNA_Tyr", "COI", "tRNA_Asp", "COII", "tRNA_Lys",
"ATP8", "ATP6", "COIII", "tRNA_Gly", "ND3", "tRNA_Arg", "ND4L",
"ND4", "tRNA_His", "ND5", "ND6", "tRNA_Glu", "CYTB", "tRNA_Thr",
"tRNA_Pro", "D_loop", "tRNA_Ser1", "tRNA_Ser2", "tRNA_Leu1",
"tRNA_Leu2"), class = "data.frame", row.names = c(NA, -1L))

setwd (tempdir())
test.out<-AnnotationBust(Accessions = "JN628859.1", Terms = mtDNAterms, Duplicates = c("tRNA_Leu","tRNA_Ser"), DuplicateInstances = c(2,2), TranslateSeqs = TRUE, TranslateCode = 2)
expect_identical(AnnotationBustOut,test.out)
})
