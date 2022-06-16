* MOST RECENT FILES/CHANGES WILL BE ON GITHUB!
global data    = "../Data/GSS"

import excel ${data}/GSS_extract.xlsx, firstrow clear

* Encode the main variables
encode satjob, gen(satiswork)     // Work satisfaction, NOT INCREASING
encode satjob1, gen(satis)        // Job satisfaction in general, INCREASING
encode toofewwk, gen(understaff)  // how often are you understaffed, NOT INCREASING
encode trynewjb, gen(tryjob)	  // How likely R make effort for new job next year
encode localnum, gen(numemp)      // Number of employees: R's work site

label list

* Replace missings
replace satiswork   = . if inlist(satiswork,1,2,3,4)
replace satis       = . if inlist(satis,1,2,3)
replace understaff  = . if inlist(understaff,1,2,3)
replace tryjob      = . if inlist(tryjob,1,2,3)
replace numemp      = . if inlist(numemp,1,2,3,4)

* Check how much data is available
tab year if !missing(understaff) & !missing(satis)
tab year if !missing(understaff) & !missing(satiswork)

*** JRG Explorations on 6/9 call

* See how many obs have different variables non-missing
tab trynewjb understaff
tab educ understaff, mi  
tab occ10, mi
tab indus10, mi

* Focus on obs which have understaffing measure
keep if !mi(understaff)

* Make string variables numeric
encode educ, gen(educ_num)   // Highest year of school completed
encode indus10, gen(indnum)  // R's industry code (NAICS 2007)
encode occ10, gen(occnum)    // R's census occupation code (2010)

* Generate binary variables
gen byte unsatisfied   = inlist(satis,4,5) // A little or very dissatisfied
gen byte verysatisfied = inlist(satis,7)   // very satisfied

destring age, replace force
gen age2 = age^2

* Generate binary variables
gen byte verylikelytryjob = tryjob == 6 if !mi(tryjob)
gen byte likelytryjob     = tryjob >= 5 if !mi(tryjob)

destring yearsjob, gen(tenure) force
replace tenure = 0.5 if yearsjob == "6-11.9 months"
replace tenure = 0 if yearsjob == "Less than 6 months"
replace tenure = . if tenure < 0

save ${data}/GSS_extract, replace

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
