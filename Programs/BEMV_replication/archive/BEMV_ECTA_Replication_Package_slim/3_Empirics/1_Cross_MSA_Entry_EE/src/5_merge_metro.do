// ***************************************************
// This script merges J2J and BDS data at the MSA level
// ---
// Alex Weinberg, November 2019
// ***************************************************


use "../output/j2j_metro", clear


* Num MSAs = 290
qui bysort msa18: gen nvals = _n == 1
count if nvals
qui drop nvals

drop if strpos(label, "Not in metropolitan area")

* Num MSAs = 256
qui bysort msa18: gen nvals = _n == 1
count if nvals
qui drop nvals


/*
!========WHICH MSAs DON'T MATCH?===========!
rename msa18 msa9
merge m:1 msa9 year using "../output/metroBDS"
qui levelsof label if _merge == 1, local(levels)
 foreach ll of local levels {
	di "`ll'"
 }
*/

// SEE CROSSWALK FOR WHY MATCH BY HAND
rename label msa_name
replace msa18 = 31100 if msa_name == "Los Angeles-Long Beach-Anaheim, CA"
replace msa18 = 42060 if msa_name == "Santa Maria-Santa Barbara, CA"
replace msa18 = 26180 if msa_name == "Urban Honolulu, HI"
replace msa18 = 14060 if msa_name == "Bloomington, IL"
replace msa18 = 36860 if msa_name == "Ottawa, IL"
replace msa18 = 29140 if msa_name == "Lafayette-West Lafayette, IN"
replace msa18 = 30500 if msa_name == "California-Lexington Park, MD"
replace msa18 = 26100 if msa_name == "Holland, MI"
replace msa18 = 19380 if msa_name == "Dayton-Kettering, OH"
replace msa18 = 44600 if msa_name == "Weirton-Steubenville, WV-OH"
replace msa18 = 37820 if msa_name == "Hermiston-Pendleton, OR"
replace msa18 = 32270 if msa_name == "Wisconsin Rapids-Marshfield, WI"
 
rename msa18 msa9
merge m:1 msa9 year using "../output/metroBDS"

drop if _merge < 3
drop _merge

// get states
merge m:1 msa9 using "../output/crosswalk_09_18_v2"
encode state, gen(STATE)
drop state
rename STATE state

* Num MSAs = 241
qui bysort msa9: gen nvals = _n == 1
count if nvals
qui drop nvals
 

rename msa9 msa
drop if missing(msa)

save "../output/metro_merged_j2j_bds", replace

