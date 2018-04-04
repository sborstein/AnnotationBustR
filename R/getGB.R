#' This is an internal function in AnnotationBustR that is amodified version of ape's read.GenBank () to get accession numbers and info from GenBank 
#' @param Accessions A vector of GenBank accession numbers.
#' @details internal function returns accessions and fasta files to measure length in AnnotationBustR
#' @return A fasta with sequence names and length info
#' @examples 
#' #a vector of 4 genbank accessions, there are two for each species.
#' genbank.accessions<-c("KP978059.1","KP978060.1","JX516105.1","JX516111.1")
#' example.getGB<-getGB(genbank.accessions)
#' @references Popescu, Andrei-Alin, Katharina T. Huber, and Emmanuel Paradis. "ape 3.0: New tools for distance-based phylogenetics and evolutionary analysis in R." Bioinformatics 28.11 (2012): 1536-1537.

getGB<-function (access.nb, seq.names = access.nb, species.names = TRUE, 
          gene.names = FALSE, as.character = FALSE) 
{
  N <- length(access.nb)
  a <- 1L
  b <- if (N > 300) 
    300L
  else N
  fl <- tempfile()
  repeat {
    URL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                  paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text")
    X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    cat(X, sep = "\n", file = fl, append = TRUE)
    if (b == N) 
      break
    a <- b + 1L
    b <- b + 300L
    if (b > N) 
      b <- N
  }
  res <- ape::read.FASTA(fl)
  if (is.null(res)) 
    return(NULL)
  attr(res, "description") <- names(res)
  if (length(access.nb) != length(res)) {
    names(res) <- gsub("\\..*$", "", names(res))
    failed <- paste(access.nb[!access.nb %in% names(res)], 
                    collapse = ", ")
    warning(paste0("cannot get the following sequences:\n", 
                   failed))
  }
  else names(res) <- access.nb
  if (as.character) 
    res <- as.character(res)
  if (species.names) {
    a <- 1L
    b <- if (N > 300) 
      300L
    else N
    sp <- character(0)
    repeat {
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                   paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
                   sep = "")
      X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE, 
                n = -1)
      sp <- c(sp, gsub(" +ORGANISM +", "", grep("ORGANISM", 
                                                X, value = TRUE)))
      if (b == N) 
        break
      a <- b + 1L
      b <- b + 300L
      if (b > N) 
        b <- N
    }
    attr(res, "species") <- gsub(" ", "_", sp)
  }
  if (gene.names) 
    warning("you used 'gene.names = TRUE': this option is obsolete; please update your code.")
  res
}
