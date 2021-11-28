## Resubmission (28/11/2021)
This is a re-submission with the following corrections, as requested:

* The length of the title has been reduced to 51 characters.

* I added more details about the package functionality and implemented methods 
in the Description text. 

* I have unwrapped all examples, and I have replaced dontrun{} with donttest{} 
in those being executable in > 5 sec.

* I have replaced options(warn=-1) with suppressWarnings().

* I have reset to user's options() (as indicated in the email) in the vignette network_description, where it was needed.  

## Resubmission (22/11/2021)
This is a re-submission with the following corrections, as requested:

* In the DESCRIPTION, the part " + file LICENCE" has been omitted and the file 
has been removed.

* Presently, there are no references about the package to add in the Description
field in the form Authors (year) <doi:10.....> or <arXiv:.....>. 

## Checking and building Windows binary package (22/11/2021)

* **R-oldrelease: R-oldrelease, currently R-4.0.5**

Status: 1 NOTE (New submission)

* **R-release: R-release, currently R-4.1.1**

Status: 1 NOTE (New submission)

* **R-devel: R-devel, to be R-4.2.0**

Status: 1 NOTE (New submission)

## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

## Resubmission (21/11/2021)
This is a re-submission. In this version, I have:

* Corrected the *'Error in xtfrm.data.frame(x) : cannot xtfrm data frames'* that
appeared in some functions, tests and vignettes. This error appeared in Debian.

## Checking and building Windows binary package

* **R-oldrelease: R-oldrelease, currently R-4.0.5**

Status: 1 NOTE (New submission)

* **R-release: R-release, currently R-4.1.1**

Status: 1 NOTE (New submission)

* **R-devel: R-devel, to be R-4.2.0**

Status: 1 NOTE (New submission)

## R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded
