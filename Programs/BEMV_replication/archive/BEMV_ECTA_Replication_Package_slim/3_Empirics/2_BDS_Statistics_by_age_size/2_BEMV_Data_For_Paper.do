clear all
set maxvar 20000
set matsize 10000
set more off
cap log close
set scheme mine
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd "$base"

local yearlb = 2013
local yearub = 2014 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CREATE NEW VARIABLES
graph drop _all

foreach xxx in "size" "age" {
// foreach xxx in "age" {

// LOAD SAMPLE
use "$working/BDS_J2J_age_size_sector", clear
keep if year >= `yearlb'
keep if year <= `yearub'

// CATEGORICAL VARIABLE
keep if type == "`xxx'"
if "`xxx'" == "age" {
	gen catvar = age_J2J 
	copydesc age_J2J catvar
	}
else if  "`xxx'" == "size" {
	gen catvar = size_J2J 
	copydesc size_J2J catvar
}

// NEW VARIABLES
/* THIS CODE IS MIRRORED IN MATLAB */
gen H 			= EEhire + NEhire 		// Hires
gen S 			= EEquit + ENsep 		// Separations
gen Poach 		= (EEhire - EEquit) 	// Poaching
gen Churn		= (H - JC) + (S - JD) 	// Churn
gen Realloc 	= (JC + JD)

gen JC_r 		= JC/N
gen JD_r 		= JD/N
gen H_r 		= H/N
gen S_r 		= S/N
gen Churn_r 	= Churn/N
gen Poach_r 	= Poach/N
gen Realloc_r 	= Realloc/N

// Print out average annual JCR and monthly hiring rate, to check these make sense: JC ~ 0.13, HR ~ 0.03
preserve
	collapse (sum) JC H N, by(year)
	gen JCR = (4*JC)/N 			// Annual job creation rate
	gen HR  = (H/3)/N 			// Monthly hiring rate
	collapse (mean) JCR HR
	list JCR HR
restore

by year: egen JCtotal = total(JC)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FIRM EXIT 
gen NX_f 		= JD_FX/N 				// Fraction of employment that exits
gen FX_f 		= FX/F 					// Fraction of firms that exit
gen JDX_f 		= JD_FX/JD 				// Fraction of job destruction due to exit
gen JC_f 		= JC/JCtotal 			// Fraction of job creation by group
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// COMPOSITION OF HIRES
gen EEhire_f 		= EEhire/H 			// Fraction of hires from employment 
gen EEquit_f 		= EEquit/S 			// Fraction of separations to employment 
gen EEhire_r 		= EEhire/N 			// EE hire rate

// COLLAPSE OVER SAMPLE: (i)Averages of rates, (ii) Add totals
collapse (mean) JC_r JD_r H_r S_r Poach_r Churn_r Realloc_r NX_f FX_f JDX_f JC_f EEhire_f EEquit_f EEhire_r ///	MEAN
		 (mean)  F E N N_Census /// 														MEAN
		 , by(catvar)
		 
// COMPUTE SHARES: Firms, Establishments, Employment
foreach X of varlist F E N N_Census {
	egen `X'_t = total(`X')
	 gen `X'_f = `X'/`X'_t
	drop `X'_t
}

keep 	catvar F_f E_f N_f N_Census_f JC_r JD_r H_r S_r Poach_r Realloc_r Churn_r NX_f FX_f JDX_f JC_f EEhire_f EEquit_f  EEhire_r
order 	catvar F_f E_f N_f N_Census_f JC_r JD_r H_r S_r Poach_r Realloc_r Churn_r NX_f FX_f JDX_f JC_f EEhire_f EEquit_f  EEhire_r

rename catvar `xxx'

gen H_r_N_f = H_r*N_f
gen H_r_N_Census_f = H_r*N_Census_f
gen S_r_N_f = S_r*N_f
gen S_r_N_Census_f = S_r*N_Census_f

export excel using "$base\output\BDS_J2J_by_`xxx'.xlsx", replace sheet("Results") firstrow(var) keepcellfmt 

// Save copy to be read into MATLAB
cd ../../2_MATLAB_code/Input_data
export excel using "BDS_J2J_by_`xxx'.xlsx", replace sheet("Results") firstrow(var) keepcellfmt 
cd "$base"

}



