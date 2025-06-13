# Submission, version 0.5.0 (13/06/2025)
* **R version 4.4.3 (2025-02-28 ucrt)**

Status: 1 NOTE
Possibly misspelled words in DESCRIPTION:
  Spineli (20:31, 21:31, 23:10, 23:59, 24:29, 25:9, 27:32, 30:10, 31:5, 33:10)
  al (21:42, 23:70, 25:20, 27:43, 31:16)
  et (21:39, 23:67, 25:17, 27:40, 31:13)
  
* **R version 4.5.0 (2025-04-11 ucrt)**

Status: 1 NOTE
Possibly misspelled words in DESCRIPTION:
  Spineli (20:31, 21:31, 23:10, 23:59, 24:29, 25:9, 27:32, 30:10, 31:5, 33:10)
  al (21:42, 23:70, 25:20, 27:43, 31:16)
  et (21:39, 23:67, 25:17, 27:40, 31:13)
* The spelling is correct!

* **R Under development (unstable) (2025-06-12 r88305 ucrt)**

Status: 1 NOTE
Possibly misspelled words in DESCRIPTION:
  Spineli (20:31, 21:31, 23:10, 23:59, 24:29, 25:9, 27:32, 30:10, 31:5, 33:10)
  al (21:42, 23:70, 25:20, 27:43, 31:16)
  et (21:39, 23:67, 25:17, 27:40, 31:13)
* The spelling is correct!

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

# Submission, version 0.4.0 (24/03/2024)
* **R version 4.2.3 (2023-03-15 ucrt)**

Status: 1 NOTE
Possibly misspelled words in DESCRIPTION:
  dendrogram (31:25)
* The spelling is correct in the indicated word!

### R CMD check results

0 errors | 0 warnings | 0 notes

R CMD check succeeded

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
