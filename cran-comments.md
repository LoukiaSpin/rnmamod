# Submission, version 0.3.0 (01/11/2022)
* **R Under development (unstable) (2022-10-11 r83083 ucrt)**

Status: OK

* **R version 4.2.2 (2022-10-31 ucrt)**

Status: OK

* **R version 4.1.3 (2022-03-10)**

Status: 1 NOTE
Possibly mis-spelled words in DESCRIPTION:
  Kullback (28:32)
  Leibler (28:41)
  heatmaps (25:23, 26:28)
  rankograms (26:67)
* The spelling is correct in all indicated words!

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

# Submission, version 0.2.0 (06/04/2022)

* **R version 4.0.5 (2021-03-31)**

Status: OK

* **R version 4.1.3 (2022-03-10)**

Status: OK

* **R version 4.2.0 alpha (2022-04-05 r82100 ucrt)**

Status: OK

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

# Resubmission (28/11/2021)
This is a re-submission with the following corrections, as requested:

* The length of the title has been reduced to 51 characters.

* I added more details about the package functionality and implemented methods 
in the Description text. 

* I have unwrapped all examples, and I have replaced dontrun{} with donttest{} 
in those being executable in > 5 sec.

* I have replaced options(warn=-1) with suppressWarnings().

* I have reset to user's options() in the vignette network_description, where it
was needed.  

### Checking and building Windows binary package 

* **R-oldrelease: R-oldrelease, currently R-4.0.5**

Status: 1 NOTE (New submission)

The following are NOT mis-spelled words in DESCRIPTION:

- Kullback (27:32)
- Leibler (27:41)
- heatmaps (24:23, 25:28)
- rankograms (25:67)

* **R-release: R-release, currently R-4.1.1**

Status: 1 NOTE (New submission)

The following are NOT mis-spelled words in DESCRIPTION:

- Kullback (27:32)
- Leibler (27:41)
- heatmaps (24:23, 25:28)
- rankograms (25:67)

* **R-devel: R-devel, to be R-4.2.0**

Status: 1 NOTE (New submission)

The following are NOT mis-spelled words in DESCRIPTION:

- Kullback (27:32)
- Leibler (27:41)
- heatmaps (24:23, 25:28)
- rankograms (25:67)

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

# Resubmission (22/11/2021)
This is a re-submission with the following corrections, as requested:

* In the DESCRIPTION, the part "+ file LICENCE" has been omitted and the file 
has been removed.

* Presently, there are no references about the package to add in the Description
field in the form Authors (year) <doi:10.....> or <arXiv:.....>. 

### Checking and building Windows binary package 

* **R-oldrelease: R-oldrelease, currently R-4.0.5**

Status: 1 NOTE (New submission)

* **R-release: R-release, currently R-4.1.1**

Status: 1 NOTE (New submission)

* **R-devel: R-devel, to be R-4.2.0**

Status: 1 NOTE (New submission)

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

# Resubmission (21/11/2021)
This is a re-submission. In this version, I have:

* Corrected the *'Error in xtfrm.data.frame(x) : cannot xtfrm data frames'* that
appeared in some functions, tests and vignettes. This error appeared in Debian.

### Checking and building Windows binary package

* **R-oldrelease: R-oldrelease, currently R-4.0.5**

Status: 1 NOTE (New submission)

* **R-release: R-release, currently R-4.1.1**

Status: 1 NOTE (New submission)

* **R-devel: R-devel, to be R-4.2.0**

Status: 1 NOTE (New submission)

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded
