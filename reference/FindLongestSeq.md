# Find the longest sequence for each species from a list of GenBank accession numbers.

Find the longest sequence for each species from a list of GenBank
accession numbers.

## Usage

``` r
FindLongestSeq(Accessions, BatchSize = 300)
```

## Arguments

- Accessions:

  A vector of GenBank accession numbers.

- BatchSize:

  Numeric. If the number of accessions is over the number provided,
  requests will be sent in batches of this amount. This is necessary for
  the NCBI servers. If you receive an HTTP 414 error, try to reduce the
  size of the batch. Default is 300.

## Value

A list of genbank accessions numbers for the longest sequence for each
taxon in a list of accession numbers.

## Details

For a set of GenBank accession numbers, this will return the longest
sequence for in the set for species.

## Examples

``` r
#a vector of 4 genbank accessions, there are two accessions for each species.
genbank.accessions<-c("KP978059.1","KP978060.1","JX516105.1","JX516111.1")
if (FALSE) { # interactive()
#returns the longest sequence respectively for the two species.
long.seq.result <- FindLongestSeq(genbank.accessions)
}
```
