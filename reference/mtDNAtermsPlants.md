# Mitochondrial DNA Search Terms for Plants

A data frame containing search terms for plant mitochondrial loci. Can
be subset for loci of interest. Columns are as follows and users should
follow the column format if they wish to add search terms using the
MergeSearchTerms function:

## Usage

``` r
mtDNAtermsPlants
```

## Format

A data frame of of 248 rows and 3 columns

- Feature: Feature name, FASTA files will be written with this name.

- Type: Type of feature, either CDS,tRNA,rRNA.

- Name: Name of synonym for a feature to search for.

## See also

[`MergeSearchTerms`](http://sborstein.github.io/AnnotationBustR/reference/MergeSearchTerms.md)
