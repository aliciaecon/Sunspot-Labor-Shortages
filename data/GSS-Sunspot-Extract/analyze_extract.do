
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

*** JRG Explorations on 6/9 call
* See how many obs have different variables non-missing
tab trynewjb understaff
tab educ understaff, mi
tab occ10, mi
tab indus10, mi

* Focus on obs which have understaffing measure
keep if !mi(understaff)

* Make string variables numeric
encode trynewjb, gen(tryjob)
replace tryjob = . if inlist(tryjob,1,2)

encode educ, gen(educ_num)
encode indus10, gen(indnum)
encode occ10, gen(occnum)

encode satjob1, gen(satis)
replace satis = . if inlist(satis,1,2)

gen byte unsatisfied = inlist(satis,3,4) if !mi(satis)
gen byte verysatisfied = inlist(satis,6) if !mi(satis)

destring age, replace force
gen age2 = age^2

gen byte verylikelytryjob = tryjob == 5 if !mi(tryjob)
gen byte likelytryjob = tryjob >= 4 if !mi(tryjob)

destring yearsjob, gen(tenure) force
replace tenure = 0.5 if yearsjob == "6-11.9 months"
replace tenure = 0 if yearsjob == "Less than 6 months"
replace tenure = . if tenure < 0

* Run regressions testing correlation of understaffing with job satisfaction and 
* probability of looking for a new job
foreach var in satis unsatisfied verysatisfied tryjob likelytryjob verylikelytryjob {
	summ `var'
	reg `var' i.understaff, vce(cluster indnum)
	reg `var' i.understaff i.educ_num, vce(cluster indnum)
	reg `var' i.understaff i.educ_num age age2, vce(cluster indnum)
	reg `var' i.understaff i.educ_num age age2 tenure, vce(cluster indnum)
	reghdfe `var' i.understaff i.educ_num age age2 tenure, vce(cluster indnum) absorb(occnum)
	reghdfe `var' i.understaff i.educ_num age age2 tenure, vce(cluster indnum) absorb(occnum indnum)
	}
