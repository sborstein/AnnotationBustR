#Find Longest Seq Test for a large number of accessions to test limit.
#Expect that the longest seq for A. christyi is accesseion
# number KT691775.1 is the longest of 1366 bp. These are non-homologous
#UCE loci, but a lot of them, so try it anyways

test_that("More than max query works to find longest seq",{
data(sysdata)
long.manual<-ape::read.GenBank(as.vector(Achrist$V1))
long.manual.res<-which(summary(long.manual)[,1]==max(as.numeric(summary(long.manual)[,1])))
expect_identical(FindLongestSeq(Achrist$V1)$Accession, names(long.manual.res))
})
