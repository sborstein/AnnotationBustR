#Test that the main function AnnotationBust is in fact working properly
#Run on one mitogenome which is missing ND5, which should return NA
#Run with duplicate genes and translated seqs
test_that("AnnotationBust Works",{
data(mtDNAterms)
data(sysdata)#Manual output of all genes. Missing ND5 gene (misc_feature). Should be missing in real run.
#test with accession=JN628859.1
setwd (tempdir())
test.out<-AnnotationBust(Accessions = "JN628859.1", Terms = mtDNAterms, Duplicates = c("tRNA_Leu","tRNA_Ser"), DuplicateInstances = c(2,2), TranslateSeqs = TRUE, TranslateCode = 2)
expect_identical(AnnotationBustOut,test.out)
})
