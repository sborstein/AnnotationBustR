#Find Longest Seq Test for a large number of accessions to test limit.
#Expect that the longest seq for A. christyi is accesseion
# number KT691775.1 is the longest of 1366 bp. These are non-homologous
#UCE loci, but a lot of them, so try it anyways

test_that("More than max query works to find longest seq",{
Achrist<-read.table("Achrist.txt")
Achrist<-Achrist$V1#1041 Accession numbers for A. burtoni
correct.seq<-ncbi_byname("Aristochromis christyi", gene="UCE")[,c(1,3:5)]#ncbi_byname from traits gets longest seq for the species for a gene
expect_identical(FindLongestSeq(Achrist)$acc_no, correct.seq$acc_no)
})

#Find Longest Seq Test for a small subset. Expect that the longest seq for A. burtoni is accesseion
# number AB915493.1 is the longest of 900 bp

test_that("Less than max query works to find longest seq",{
  Aburt<-read.table("Aburt.txt")
  Aburt<-Aburt$V1#8 Accession numbers for A. burtoni
  correct.seq<-ncbi_byname("Haplochromis burtoni", gene="cytb")[,c(1,3:5)]#ncbi_byname from traits gets longest seq for the species for a gene
  expect_identical(FindLongestSeq(Aburt)$acc_no, correct.seq$acc_no)
})