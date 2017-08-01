#' Mitochondrial DNA Search Terms for Animals
#'
#' A data frame containing search terms for animal mitochondrial loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 174 rows and 3 columns
#' \itemize{
#'   \item Locus: Locus name, FASTA files will be written with this name
#'   \item Type: Type of subsequence, either CDS,tRNA,rRNA,misc_RNA, or D-loop
#'   \item Name:Name of synonym for a locus to search for
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"mtDNAterms"

#' Chloroplast DNA (cpDNA) Search Terms
#'
#' A data frame containing search terms for Chloroplast loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 201 rows and 3 columns
#' \itemize{
#'   \item Locus: Locus name, FASTA files will be written with this name
#'   \item Type: Type of subsequence, either CDS,tRNA,rRNA, or misc_RNA
#'   \item Name:Name of synonym for a locus to search for
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"cpDNAterms"

#' Ribosomal DNA (rDNA) Search Terms
#'
#' A data frame containing search terms for ribosomal RNA loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 7 rows and 3 columns
#' \itemize{
#'   \item Locus: Locus name, FASTA files will be written with this name
#'   \item Type: Type of subsequence, either rRNA or misc_RNA
#'   \item Name:Name of synonym for a locus to search for
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"rDNAterms"