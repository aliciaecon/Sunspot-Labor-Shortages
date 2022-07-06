// ***************************************************
// This script cleans BDS data at the MSA level
// ---
// Alex Weinberg, November 2019
// ***************************************************

import delimited "../input/metroBDS.csv", clear
rename year2 year 


// PACKAGE FOR MERGE WITH J2J --------------------------------------------------
drop if year < 2000 // j2j data starts in 2000

rename msa msa9

// SAVE 
sort msa9 year
save "../output/metroBDS", replace
