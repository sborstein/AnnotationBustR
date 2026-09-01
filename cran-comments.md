## Resubmission/Update
This is a resubmission of a package update of AnnotationBustR from v1.3 to v2.0. In this version I have:

* Update package dependencies to allow access to more current data.
* Added proper Authors@R field in DESCRIPTION.
* Removed calls to structure() using deprecated special names in unit tests.
* Note some function have unit tests using testthat for functions in which the examples are wrapped with /donttest. It is not feasible to unwrap these functions from /donttest as the functions and examples connect to a remote database and take longer than 5 seconds to run. The function AnnotationBust and FindLongestSeq have testthat test file with similar names. I have checked the examples of these functions that are wrapped in /donttest using R CMD CHECK --run-donttest, and all ran without error.
* win-builder (oldrelease) checks suggest some URLs may be invalid. We have checked the URLs and they are correct and can be reached via web browser and are valid.


## Test environments
* local Windows 11 install, R 4.6.1
* Ubuntu 24.04.3 LTS on GitHub Actions(devel, release, oldrelease)
* win-builder (devel, release, and oldrelease)
* macOS 15 on GitHub Actions (release)

## R CMD check results
* There were no ERRORS or WARNINGS