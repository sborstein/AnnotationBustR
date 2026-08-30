# Ribosomal DNA (rDNA) Search Terms

A data frame containing search terms for ribosomal RNA loci. Can be
subset for loci of interest. Columns are as follows and users should
follow the column format if they wish to add search terms using the
MergeSearchTerms function:

## Usage

``` r
rDNAterms
```

## Format

A data frame of of 14 rows and 3 columns

- Feature: Feature name, FASTA files will be written with this name.

- Type: Type of feature, either rRNA or misc_RNA.

- Name: Name of synonym for a feature to search for.

## See also

[`MergeSearchTerms`](http://sborstein.github.io/AnnotationBustR/reference/MergeSearchTerms.md)
