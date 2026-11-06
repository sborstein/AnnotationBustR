#Test to see if merge search terms works properly.
#Will merge rDNA and Mitogenome terms together, which will be included in the package as data. 
#Do one to test if it merges and another to see if it sorts
test_that("merging search terms works", {
  data("mtDNAterms")#load mito search terms
  data("rDNAterms")#load rDNA search terms
  #No Sorting
  reg_merged<-rbind(mtDNAterms, rDNAterms)#simple rbind
  reg.test.merge<-MergeSearchTerms(mtDNAterms, rDNAterms, SortGenes = FALSE)#run the function w/out sorting
  expect_identical(dim(reg.test.merge), dim(reg_merged))#test to see if same dimensions
  expect_identical(reg.test.merge[1,], reg_merged[1,])#test to see the first row is the same 
  expect_identical(reg.test.merge[dim(reg.test.merge)[[1]],], reg_merged[dim(reg_merged)[[1]],])#test to see the last row is the same 
  #With Sorting
  sorted_merged<-reg_merged[order(reg_merged$Feature),]#order the manual ones
  row.names(sorted_merged)<-1:dim(sorted_merged)[[1]]#re-number
  sort.test.merge<-MergeSearchTerms(mtDNAterms, rDNAterms, SortGenes = TRUE)#run the function with sorting
  expect_identical(dim(sort.test.merge), dim(sorted_merged))#test to see if same dimensions
  expect_identical(sort.test.merge[1,], sorted_merged[1,])#test to see the first row is the same 
  expect_identical(sort.test.merge[dim(sort.test.merge)[[1]],], sorted_merged[dim(sorted_merged)[[1]],])#test to see the last row is the same
  
  #Check Fail Conditions
  #Different dimension terms
  expect_error(MergeSearchTerms(mtDNAterms, rDNAterms[,-3], SortGenes = FALSE))
  #Wrong format
  expect_error(MergeSearchTerms(mtDNAterms, as.matrix(rDNAterms), SortGenes = FALSE))
  #Deviations from 3 or 4 columns
  expect_error(MergeSearchTerms(mtDNAterms[,1:2], rDNAterms[,1:2], SortGenes = FALSE))
  #Deviations in Column Names
  mtDNAterms.colChange <- mtDNAterms#store as new object to change column names
  rDNAterms.colChange <- rDNAterms#store as new object to change column  names
  colnames(mtDNAterms.colChange) <-c("Locus","Type","Name")
  colnames(rDNAterms.colChange) <-c("Locus","Type","Name")
  expect_error(MergeSearchTerms(mtDNAterms.colChange, rDNAterms.colChange, SortGenes = FALSE))
}
)

