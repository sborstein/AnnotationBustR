## Resubmission
This is a resubmission of a package update of AnnotationBustR from v1.2 to v1.3. In this version I have:

*Addressed issues where unit tests failed on CRAN due to not being able to successfully connect to a remote database. These tests will now skip if a connection can not be made. Please note that these tests ran successfully on Travis-CI as well as R WinBuilder and on my local Windows 10 and Mac OS environments, but previous submissions have had issues connecting to the database during CRAN checks. Note some examples, such as those in the AnnotationBust and FindLongestSeq functions, connect to remote databases.
*Note some function have unit tests using testthat for functions in which the examples are wrapped with /donttest. It is not feasible to unwrap these functions from /donttest as the functions and examples connect to a remote database and take longer than 5 seconds to run. The function AnnotationBust and FindLongestSeq have testthat test file with similar names. I have checked the examples of these functions that are wrapped in /donttest using R CMD CHECK --run-donttest, and all ran without error.

## Test environments
* local OS X install, R 4.0.2
* ubuntu 12.04 (on travis-ci), R 4.0.2
* local Windows 10 install, R 4.0.2
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes
For R CMD check ran with ubuntu, OS X, win-builder (dev), and win-builder (release)
