// ***************************************************
// This script make a crosswalk for the MSA codes
// J2J uses the 2018 MSA codes: https://lehd.ces.census.gov/data/schema/j2j_latest/lehd_public_use_schema.html#geo_level
// BDS uses the 2009 MSA codes: https://www.census.gov/content/dam/Census/programs-surveys/business-dynamics-statistics/BDS_Codebook.pdf
//
// This script merges creates a crosswalk so that we can merge J2J and BDS. 
// ---
// Alex Weinberg, November 2019
// ***************************************************

// 2018 MSA CODES ----------------------------------------------

import excel "../input/msa18.xls", clear
drop in 1/3
drop in 1916/l

rename E size
rename D msa_name
rename I state
destring A, gen(msa18)

// keep only one obs per MSA
bysort msa18: gen nvals = _n == 1
di "Number of unique CBSA regions"
count if nvals
keep if nvals

keep msa18 msa_name state size
save "../output/msa18_codes", replace

// 2009 MSA CODES ----------------------------------------------

import excel "../input/msa9.xls", clear

drop in 1/4
drop in 1863/l

rename E size
rename D msa_name
rename J state
destring A, gen(msa9)

// keep only one obs per MSA
bysort msa9: gen nvals = _n == 1
di "Number of unique CBSA regions"
count if nvals
keep if nvals

keep msa9 msa_name state size
save "../output/msa9_codes", replace

// MERGE MSA CODES ---------------------------------------------

merge 1:1 msa_name state using "../output/msa18_codes"

// All the successful merges have the same MSA code 2009/2018
count if _merge==3
count if _merge==3 & msa9 == msa18


// Same MSA code, different name of the MSA region
qui levelsof msa9 if _merge < 3, local(levels)
 foreach ll of local levels {
	qui count if msa18 == `ll'
	if (r(N) == 1) {
		list msa9 msa_name msa18 if msa18 == `ll' | msa9 == `ll'
		qui replace msa18 = msa9 if msa9 == `ll' 
		qui drop if msa18 == `ll'  & missing(msa9)
		qui replace _merge = 3 if msa9 == `ll' 
	}	
 }

 // MATCH REMAINING CODES BY HAND (different name+ different code) -------------
/* One concern is that these MSAs are changing that that drives the changes in J2J and BDS 
e.g. new_firms jumps because the MSA now includes a new city, not because acutally jumps
I'm not so worried because 2006-2009 so short period and I think difference is naming convention, not actual different area
I will look into later.
*/
sort state msa_name
count if missing(msa9) | missing(msa18)
list msa9 msa_name if !missing(msa9) & missing(msa18)

replace msa18 = msa9 if msa_name == "Prescott, AZ"  // this is what merges
drop if msa_name == "Prescott Valley-Prescott, AZ"

replace msa18 = 31080 if msa_name == "Los Angeles-Long Beach-Santa Ana, CA"
drop if msa_name == "Los Angeles-Long Beach-Anaheim, CA"

replace msa18 = 42200 if msa_name == "Santa Barbara-Santa Maria-Goleta, CA"
drop if msa_name == "Santa Maria-Santa Barbara, CA"

replace msa18 = 46520 if msa_name == "Honolulu, HI"
drop if msa_name == "Urban Honolulu, HI"

replace msa18 = 14010 if msa_name == "Bloomington-Normal, IL"
drop if msa_name == "Bloomington, IL"

replace msa18 = 36837 if msa_name == "Ottawa-Streator, IL"
drop if msa_name == "Ottawa, IL"

replace msa18 = 29200 if msa_name == "Lafayette, IN"
drop if msa_name == "Lafayette-West Lafayette, IN"

replace msa18 = 15680 if msa_name == "Lexington Park, MD"
drop if msa_name == "California-Lexington Park, MD"

replace msa18 = 26090 if msa_name == "Holland-Grand Haven, MI"
drop if msa_name == "Holland, MI"

replace msa18 = 19430 if msa_name == "Dayton, OH"
drop if msa_name == "Dayton-Kettering, OH"

replace msa18 = 48260 if msa_name == "Steubenville-Weirton, OH-WV"
drop if msa_name == "Weirton-Steubenville, WV-OH"

replace msa18 = 25840 if msa_name == "Pendleton-Hermiston, OR"
drop if msa_name == "Hermiston-Pendleton, OR"

replace msa18 = 49220 if msa_name == "Marshfield-Wisconsin Rapids, WI"
drop if msa_name == "Wisconsin Rapids-Marshfield, WI"

drop _merge
di "Number of MSAs"
count if !missing(msa18) & !missing(msa9) & size=="Metropolitan Statistical Area"
di "Number of microSAs"
count if !missing(msa18) & !missing(msa9) & size!="Metropolitan Statistical Area"

keep  if !missing(msa18) & !missing(msa9) 
save "../output/crosswalk_09_18_v2", replace
