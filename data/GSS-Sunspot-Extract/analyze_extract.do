
import excel GSS.xlsx, firstrow clear

* Group the relevant variables
egen understaff   = group(toofewwk)
egen jobsatis     = group(satjob)
egen jobsatis2    = group(satjob1)


*net install gr0034.pkg
* https://stats.oarc.ucla.edu/stata/faq/how-do-i-assign-the-values-of-one-variable-as-the-value-labels-for-another-variable/

labmask jobsatis, values(satjob)
labmask jobsatis2, values(satjob1)
labmask understaff, values(toofewwk)

label list

* Keep data for which we have understaffed and job/work satisfaction 
replace jobsatis    = . if inlist(jobsatis, 1,2,3,4)
replace understaff  = . if inlist(understaff, 1,2,3)
replace jobsatis2   = . if inlist(jobsatis2, 1,2,3)

* Check how much data is available
count if !missing(jobsatis)
count if !missing(jobsatis2)
count if !missing(understaff)

tab year if !missing(understaff) & !missing(jobsatis2)
tab year if !missing(understaff) & !missing(jobsatis)

* Check the correlations
pwcorr jobsatis understaff
pwcorr jobsatis2 understaff
pwcorr jobsatis*
