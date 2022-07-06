clear all
set maxvar 20000
set matsize 10000
set more off
cap log close
set scheme mine
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd "$base"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CONTENTS:
/* BDS */
// IMPORT BDS DATA - AGE
// IMPORT BDS DATA - SIZE
// IMPORT BDS DATA - SECTOR
/* J2J */
// IMPORT J2J DATA - AGE
// IMPORT J2J DATA - SIZE
// IMPORT J2J DATA - SECTOR
/* Combine */
// COMBINE - AGE
// COMBINE - SIZE
// COMBINE - SECTOR
// APPEND - AGE+SIZE+SECTOR
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT BDS DATA - AGE
import delimited using "$base/raw_data_BDS/bds_f_age_release.csv",  clear
rename year2 year
// BDS CATEGORIES
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gen 	age = .
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace age = 1  if fage4 == "a) 0"  			// J2J: {0-1}
replace age = 2  if fage4 == "b) 1"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace age = 3  if fage4 == "c) 2"	 			// J2J: {2-3}
replace age = 4  if fage4 == "d) 3"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace age = 5  if fage4 == "e) 4"	 			// J2J: {4-5}
replace age = 6  if fage4 == "f) 5"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace age = 7  if fage4 == "g) 6 to 10"	 	// J2J: {6-10}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace age = 8  if fage4 == "h) 11 to 15"	 	// J2J: {11+}
replace age = 9  if fage4 == "i) 16 to 20"	
replace age = 10 if fage4 == "j) 21 to 25"	
replace age = 11 if fage4 == "k) 26+"	
replace age = 12 if fage4 == "l) Left Censored"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Make conformable with J2J
recode age (1/2 = 1) (3/4 = 2) (5/6 = 3) (7 = 4) (8/12 = 5), gen(age_J2J)
label define label_firmage 1 "0-1 Years" 2 "2-3 Years" 3 "4-5 Years" 4 "6-10 Years" 5 "11+ Years"
label values age_J2J label_firmage

// Create variables consistently
gen JC		= job_creation
gen JD 		= job_destruction
gen JD_EX 	= job_destruction_deaths 		// Establishment exit
gen JD_FX 	= firmdeath_emp 				// Firm exit
gen N  		= emp
gen F  		= firms
gen E  		= estabs
gen EX 		= estabs_exit
gen FX  	= firmdeath_firms

// COLLAPSE ACROSS AGE GROUPS 
sort year age_J2J
collapse (sum) JC JD JD_EX JD_FX N E EX F FX, by(year age_J2J)

gen sec_BDS		= .
gen size_J2J 	= .
gen type 		= "age" 

keep 	year JC JD JD_EX JD_FX N E EX F FX sec_BDS age_J2J size_J2J type
order 	year JC JD JD_EX JD_FX N E EX F FX sec_BDS age_J2J size_J2J type

// Make BDS flow measures quarterly
foreach var of varlist JC JD JD_EX JD_FX JD_EX EX FX {
	replace `var' = `var'/4
}

// Annual, 1977 - 2014
format year %ty
sort year age_J2J
keep if year >= 2000 & year <= 2016 
save "$working/BDS_by_age", replace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT BDS DATA - SIZE
import delimited using "$base/raw_data_BDS/bds_f_isz_release.csv",  clear

rename year2 year

// BDS CATEGORIES
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gen 	size = .
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace size = 1  if ifsize == "a) 1 to 4"	 		// J2J: {0-19}
replace size = 2  if ifsize == "b) 5 to 9"	
replace size = 3  if ifsize == "c) 10 to 19"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace size = 4  if ifsize == "d) 20 to 49"		// J2J: {20-49}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace size = 5  if ifsize == "e) 50 to 99"		// J2J: {50-249}
replace size = 6  if ifsize == "f) 100 to 249"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace size = 7  if ifsize == "g) 250 to 499"		// J2J: {250-499}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replace size = 8  if ifsize == "h) 500 to 999"	 	// J2J: {500+}
replace size = 9  if ifsize == "i) 1000 to 2499"	
replace size = 10 if ifsize == "j) 2500 to 4999"	
replace size = 11 if ifsize == "k) 5000 to 9999"	
replace size = 12 if ifsize == "l) 10000+"	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Make conformable with J2J
recode size  (1/3 = 1) (4 = 2) (5/6 = 3) (7 = 4) (8/12 = 5), gen(size_J2J)
label define label_firmsize	1 "0-19 Employees" 2 "20-49 Employees" 3 "50-249 Employees" 4 "250-499 Employees" 5 "500+ Employees"
label values size_J2J label_firmsize

// Create variables consistently
gen JC		= job_creation
gen JD 		= job_destruction
gen JD_EX 	= job_destruction_deaths 		// Establishment exit
gen JD_FX 	= firmdeath_emp 				// Firm exit
gen N  		= emp
gen F  		= firms
gen E  		= estabs
gen EX 		= estabs_exit
gen FX  	= firmdeath_firms

// COLLAPSE ACROSS SIZE GROUPS 
sort year size_J2J
collapse (sum) JC JD JD_EX JD_FX N E EX F FX, by(year size_J2J)

replace JD_FX = JD_EX 	if (JD_FX==0)
replace FX 	  = EX	 	if (FX==0)

gen sec_BDS		= .
gen age_J2J 	= .
gen type 		= "size" 

keep 	year JC JD JD_EX JD_FX N  E EX FX F sec_BDS age_J2J size_J2J type
order 	year JC JD JD_EX JD_FX N  E EX FX F sec_BDS age_J2J size_J2J type

// Make BDS flow measures quarterly
foreach var of varlist JC JD JD_EX JD_FX JD_EX EX FX {
	replace `var' = `var'/4
}
// Annual, 1977 - 2016
format year %ty
sort year size_J2J
keep if year >= 2000 & year <= 2016
save "$working/BDS_by_size", replace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT BDS DATA - SECTOR
import delimited using "$base/raw_data_BDS/bds_f_sic_release.csv",  clear

rename year2 year

// See: Bottom of page here -> https://www.census.gov/ces/dataproducts/bds/definitions.html#sector
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gen 	sec_BDS = .
replace sec_BDS = 1 if 	sic1 == 7 								//  "Agriculture"
replace sec_BDS = 2 if 	sic1 == 10 								//  "Mining"
replace sec_BDS = 3 if 	sic1 == 15 								//  "Construction"
replace sec_BDS = 4 if 	sic1 == 20 								//  "Manufacturing"
replace sec_BDS = 5 if 	sic1 == 40 								//  "Transportation & Utilities"
replace sec_BDS = 6 if 	sic1 == 50 								//  "Wholesale Trade"
replace sec_BDS = 7 if 	sic1 == 52 								//  "Retail Trade"
replace sec_BDS = 8 if 	sic1 == 60 								//  "Finance, Insurance, Real Estate"
replace sec_BDS = 9 if 	sic1 == 70 								//  "Services"

label define label_sec 1 "Agriculture" 2 "Mining" 3 "Construction" 4 "Manufacturing" 5 "Transportation & Utilities" ///
					   6 "Wholesale Trade" 7 "Retail Trade" 8 "FIRE" 9 "Services"
label values sec_BDS label_sec

// Create variables consistently
gen JC		= job_creation
gen JD 		= job_destruction
gen JD_EX 	= job_destruction_deaths
gen JD_FX 	= firmdeath_emp
gen N  		= emp
gen F  		= firms
gen E  		= estabs
gen EX 		= estabs_exit
gen FX  	= firmdeath_firms

sort year sec_BDS
collapse (sum) JC JD JD_EX JD_FX N E EX F FX, by(year sec_BDS)

gen age_J2J 	= .
gen size_J2J 	= .
gen type 		= "sector" 


keep 	year JC JD JD_EX JD_FX N  E EX F FX sec_BDS age_J2J size_J2J type 
order 	year JC JD JD_EX JD_FX N  E EX F FX sec_BDS age_J2J size_J2J type 

// Make BDS flow measures quarterly
foreach var of varlist JC JD JD_EX JD_FX EX FX {
	replace `var' = `var'/4
}

// Annual, 1977 - 2014
format year %ty
sort year sec
keep if year >= 2000 & year <= 2014 
save "$working/BDS_by_sector", replace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT J2J DATA - AGE
// xlsx file created by hand from download: 
	// j2j_us_all_aggregations_s, j2j3_us_all_aggregations_s
import excel using "$base/raw_data_J2J/j2j_AGE.xlsx", sheet("for_Stata_level") firstrow clear

drop if firmage == "N"
drop if firmage == "0"
destring firmage, force replace
label define label_firmage 1 "0-1 Years" 2 "2-3 Years" 3 "4-5 Years" 4 "6-10 Years" 5 "11+ Years"
label values firmage label_firmage
rename firmage age_J2J 

drop EEHire EESep NEHire ENSep J2JHire J2JSep
rename MainB 		N_Census
rename EEHireBEMV 	EEhire
rename EESepBEMV	EEquit
rename NEHireBEMV	NEhire 	
rename ENSepBEMV 	ENsep 	
destring NEhire ENsep, force replace

// Take mean within each age cell and year, so averaging across 4 quarters
collapse (mean) EEhire EEquit NEhire ENsep N_Census, by(year age_J2J)
// Quarterly measures

gen sec_BDS		= .
gen size_J2J 	= .
gen type 		= "age" 

keep 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type
order 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type

format year %ty
sort year age_J2J
keep if year >= 2000 & year <= 2014 
save "$working/J2J_by_age", replace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT J2J DATA - SIZE
// xlsx file created by hand from download: 
	// j2j_us_all_aggregations_s, j2j3_us_all_aggregations_s
import excel using "$base/raw_data_J2J/j2j_SIZE.xlsx", sheet("for_Stata_level") firstrow clear

drop if firmsize == "N"
drop if firmsize == "0"
destring firmsize, force replace
label define label_firmsize	1 "0-19 Employees" 2 "20-49 Employees" 3 "50-249 Employees" 4 "250-499 Employees" 5 "500+ Employees"
label values firmsize label_firmsize
rename firmsize size_J2J

drop EEHire EESep NEHire ENSep J2JHire J2JSep
rename MainB 		N_Census
rename EEHireBEMV 	EEhire
rename EESepBEMV	EEquit
rename NEHireBEMV	NEhire 	
rename ENSepBEMV 	ENsep 	
destring NEhire ENsep, force replace

collapse (mean) EEhire EEquit NEhire ENsep (mean) N_Census, by(year size_J2J)

gen sec_BDS		= .
gen age_J2J 	= .
gen type 		= "size" 

keep 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type
order 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type

format year %ty
sort year size_J2J
keep if year >= 2000 & year <= 2014 
save "$working/J2J_by_size", replace

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// IMPORT J2J DATA - SEC
// xlsx file created by hand from download: 
	// j2j_us_all_aggregations_s, j2j3_us_all_aggregations_s
import excel using "$base/raw_data_J2J/j2j_SECTOR.xlsx", sheet("for_Stata_level") firstrow clear

drop if industry == "N"
drop if industry == "00" 		// ALL NAICS SECTORS

gen sec_BDS = .
replace sec_BDS = 1 if industry == "11"
replace sec_BDS = 2 if industry == "21"
replace sec_BDS = 5 if industry == "22"
replace sec_BDS = 3 if industry == "23"
replace sec_BDS = 4 if industry == "31-33"
replace sec_BDS = 6 if industry == "42"
replace sec_BDS = 7 if industry == "44-45"
replace sec_BDS = 5 if industry == "48-49"
replace sec_BDS = 8 if industry == "51"
replace sec_BDS = 8 if industry == "52"
replace sec_BDS = 8 if industry == "53"
replace sec_BDS = 9 if industry == "54"
replace sec_BDS = 9 if industry == "55"
replace sec_BDS = 9 if industry == "56"
replace sec_BDS = 9 if industry == "61"
replace sec_BDS = 9 if industry == "62"
replace sec_BDS = 9 if industry == "71"
replace sec_BDS = 9 if industry == "72"
replace sec_BDS = 9 if industry == "81"
drop if industry == "92"

label define label_sec 1 "Agriculture" 2 "Mining" 3 "Construction" 4 "Manufacturing" 5 "Transportation & Utilities" ///
					   6 "Wholesale Trade" 7 "Retail Trade" 8 "FIRE" 9 "Services"
label values sec_BDS label_sec

// J2J USES NAICS2:
// https://www.census.gov/cgi-bin/sssd/naics/naicsrch?chart=2017
//  Sector	Description
// 11	Agriculture, Forestry, Fishing and Hunting 			// AGRICULTURE 	(sec=1)
// 21	Mining, Quarrying, and Oil and Gas Extraction 		// MINING  		(sec=2)
// 22	Utilities 											//  ... (sec=5)
// 23	Construction 										// CONSTRUCTION (sec=3)
// 31-33	Manufacturing 									// MANUF 		(sec=4)
// 42	Wholesale Trade										// W-TRADE 		(sec=6)
// 44-45	Retail Trade 									// R-TRADE 		(sec=7)
// 48-49	Transportation and Warehousing 					// ... (sec=5)
// 51	Information 										// .....(sec=8)
// 52	Finance and Insurance 								// .....(sec=8)
// 53	Real Estate and Rental and Leasing 					// .....(sec=8)
// 54	Professional, Scientific, and Technical Services 	// .........(sec=9)
// 55	Management of Companies and Enterprises 			// .........(sec=9)
// 56	Administrative and Support and Waste Management and Remediation Services
// 61	Educational Services 								// .........(sec=9)
// 62	Health Care and Social Assistance 					// .........(sec=9)
// 71	Arts, Entertainment, and Recreation 				// .........(sec=9)
// 72	Accommodation and Food Services 					// .........(sec=9)
// 81	Other Services (except Public Administration) 		// .........(sec=9)
// 92	Public Administration 	// Drop

drop EEHire EESep NEHire ENSep J2JHire J2JSep
rename MainB 		N_Census
rename EEHireBEMV 	EEhire
rename EESepBEMV	EEquit
rename NEHireBEMV	NEhire 	
rename ENSepBEMV 	ENsep 	
destring NEhire ENsep, force replace

collapse (sum) EEhire EEquit NEhire ENsep (mean) N_Census, by(year sec_BDS)

gen size_J2J	= .
gen age_J2J 	= .
gen type 		= "sec" 

keep 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type
order 	year N_Census EEhire EEquit NEhire ENsep sec_BDS age_J2J size_J2J type

format year %ty
sort year sec_BDS
keep if year >= 2000 & year <= 2014 
save "$working/J2J_by_sector", replace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// COMBINE - AGE - BDS + J2J
use "$working/BDS_by_age", clear
merge 1:1 year age_J2J using "$working/J2J_by_age"
save "$working/BDS_J2J_by_age", replace

// COMBINE - SIZE - BDS + J2J
use "$working/BDS_by_size", clear
merge 1:1 year size_J2J using "$working/J2J_by_size"
save "$working/BDS_J2J_by_size", replace

// COMBINE - SECTOR - BDS + J2J + Kehrig
use "$working/BDS_by_sector", clear
merge 1:1 year sec_BDS using "$working/J2J_by_sector"
save "$working/BDS_J2J_by_sector", replace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// APPEND - AGE+SIZE+SECTOR
use "$working/BDS_J2J_by_size", clear
append using "$working/BDS_J2J_by_age"
append using "$working/BDS_J2J_by_sector"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* LABEL AGAIN */
label drop _all
label define label_sec 1 "Agriculture" 2 "Mining" 3 "Construction" 4 "Manufacturing" 5 "Transportation & Utilities" ///
					   6 "Wholesale Trade" 7 "Retail Trade" 8 "FIRE" 9 "Services"
label values sec_BDS label_sec

label define label_firmage 1 "0-1 Years" 2 "2-3 Years" 3 "4-5 Years" 4 "6-10 Years" 5 "11+ Years"
label values age_J2J label_firmage

label define label_firmsize	1 "0-19 Employees" 2 "20-49 Employees" 3 "50-249 Employees" 4 "250-499 Employees" 5 "500+ Employees"
label values size_J2J label_firmsize

// J2J - 2000 - 2018
// BDS - 1974 - 2014 	(can extend to 2016 using more data on website)
keep if year >= 2000 & year <= 2016
drop _merge

order year type sec_BDS age_J2J size_J2J N E EX FX F JC JD JD_FX JD_EX N_Census EEhire EEquit NEhire ENsep 
sort year  type sec_BDS age_J2J size_J2J

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
save "$working/BDS_J2J_age_size_sector", replace
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

