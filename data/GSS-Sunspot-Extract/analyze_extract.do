import excel GSS.xlsx, firstrow clear

* Group the relevant variables
encode toofewwk, gen(understaff)
encode satjob, gen(jobsat)
encode satjob1, gen(jobsat2)
encode trynewjb, gen(tryjob)

replace jobsat      = . if inlist(jobsat, 1,2,3,4)
replace jobsat2     = . if inlist(jobsat2, 1,2,3)
replace understaff  = . if inlist(understaff, 1,2,3)
replace tryjob      = . if inlist(tryjob, 1,2,3)

* Check how much data is available
count if !missing(jobsat)
count if !missing(jobsat2)
count if !missing(understaff)
count if !missing(tryjob)

tab year if !missing(understaff) & !missing(jobsat2)
tab year if !missing(understaff) & !missing(jobsat)
tab year if !missing(understaff) & !missing(tryjob)

* Reorder understaff, so increasing in how often workplace is understaffed
gen understaff2     = 1 if understaff == 4
replace understaff2 = 2 if understaff == 6
replace understaff2 = 3 if understaff == 7
replace understaff2 = 4 if understaff == 5
drop understaff
rename understaff2 understaff

* Check the correlations
pwcorr jobsat* understaff tryjob

* Drop 89 or older + missing
destring age, force replace
