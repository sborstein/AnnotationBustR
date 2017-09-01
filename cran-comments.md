## Resubmission
This is a resubmission of a package update of AnnotationBustR from v1.0 to v1.1. In this version I have:

* Added a URL for GenBank
* Added a URL rather than a DOI linking to a description of FASTA. FASTA format originated from a now defunct software program which does not mention the format by name in the paper describing the program, yet the format has survived and is still widely used to date.

## Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.1
* win-builder (devel and release), R 3.4.1

## R CMD check results
0 errors | 0 warnings | 1 notes
For R CMD check ran with ubuntu, OS X,win-builder (dev), and win-builder (release)

Found the following (possibly) invalid URLs: URL: http://cran.rstudio.com/web/packages/AnnotationBustR/index.html https://cran.rstudio.com/web/packages/AnnotationBustR/index.html

These are false positives as they are not URLs to a package but rather to markdown badges for build version and download stats.
