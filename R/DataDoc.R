#' Mitochondrial DNA Search Terms for Animals
#'
#' A data frame containing search terms for animal mitochondrial loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 254 rows and 3 columns
#' \itemize{
#'   \item Feature: Feature name, FASTA files will be written with this name.
#'   \item Type: Type of feature, CDS,tRNA,rRNA,misc_feature, or D-loop.
#'   \item Name: Name of synonym for a feature to search for.
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"mtDNAterms"

#' Mitochondrial DNA Search Terms for Plants
#'
#' A data frame containing search terms for plant mitochondrial loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 248 rows and 3 columns
#' \itemize{
#'   \item Feature: Feature name, FASTA files will be written with this name.
#'   \item Type: Type of feature, either CDS,tRNA,rRNA.
#'   \item Name: Name of synonym for a feature to search for.
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"mtDNAtermsPlants"

#' Chloroplast DNA (cpDNA) Search Terms
#'
#' A data frame containing search terms for Chloroplast loci. Can be subset for loci of 
#'  interest. Columns are as follows and users should follow the column format if they wish to
#'  add search terms using the MergeSearchTerms function:
#'
#' @format A data frame of of 364 rows and 3 columns
#' \itemize{
#'   \item Feature: Feature name, FASTA files will be written with this name.
#'   \item Type: Type of feature, either CDS,tRNA,rRNA.
#'   \item Name: Name of synonym for a feature to search for.
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
#'   \item Feature: Feature name, FASTA files will be written with this name.
#'   \item Type: Type of feature, either rRNA or misc_RNA.
#'   \item Name: Name of synonym for a feature to search for.
#' }
#' 
#' @seealso \code{\link{MergeSearchTerms}}
"rDNAterms"