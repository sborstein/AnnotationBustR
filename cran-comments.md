## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (devel and release) , R 3.3.1

## R CMD check results
0 errors | 0 warnings |  notes
For R CMD check ran with ubuntu, OS X, and win-builder (release)

0 errors | 0 warnings |  2 notes
For R CMD check ran with win-builder (devel)

* running examples for arch 'i386' ... [30s] NOTE
Examples with CPU or elapsed time > 10s
               user system elapsed
AnnotationBust 0.55   0.07   28.36

* running examples for arch 'x64' ... [28s] NOTE
Examples with CPU or elapsed time > 10s
               user system elapsed
AnnotationBust  0.5   0.16   26.89

The function AnnotationBust and its examples connect with the ACNUC database from the seqinr package remotely and can take between 20-35s to run.

