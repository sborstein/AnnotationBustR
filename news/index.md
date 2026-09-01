# Changelog

## AnnotationBustR 2.0

### NEW FEATURES

- AnnotationBustR 2.0 now uses rentrez to provide access to sequences.
  This offers access to the latest versions of GenBank and connection
  stability.
- AnnotationBustR 2.0 now uses various BioConductor packages as
  dependencies for parsing annotations, providing speed increases and
  access to more feature types.
- Users no longer have to provide a translation code when requesting
  coding sequences be translated.
- Users can now request coding sequences to be returned as DNA or the
  translated sequence or both.
- Accession tables written by AnnotationBustR now contain a citation for
  the generation of sequence data, allowing users to properly cite and
  give credit to the researchers who generated the sequence.

### UPDATED FEATURES

- Updated vignette
- Updated unit tests to behave better if a database connection cannot be
  made
