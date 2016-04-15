#Test to see if annotation positions are found correctly

test_that("Getting Sequence Positions from Annotations works",{
numbs<-c("FJ706343","FJ706292")
rdna.test.positions<-GetSeqPos(numbs, rDNA.Genes, duplicate.genes = NULL)
rdna.manual.positions<-read.csv("rdna_manualPositions.csv", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
class(rdna.manual.positions)<-append(class(rdna.manual.positions),"Annot.Pos")#give the manual positions class Annot.Pos to check identicality
expect_that(rdna.test.positions, is_a("Annot.Pos"))#Check that it is of proper class Annot.Pos
expect_identical(dim(rdna.test.positions), dim(rdna.manual.positions))#test to see if same dimensions
expect_identical(rdna.test.positions[,1:2], rdna.manual.positions[,1:2])#test to see if they are the same
expect_identical(as.numeric(unlist(rdna.test.positions[,3:12])), as.numeric(unlist(rdna.manual.positions[,3:12])))#test to see if they are the same
})