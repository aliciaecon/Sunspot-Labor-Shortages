// ***************************************************
// 1. Compute j2j hire rate and firm entry rate for each metro.
// 2. Compute ratio btwn the above rates in 2006:Q2 and 2009:Q2 (follows HHKM 2018 AEJ Macro)
// 3. Plot the correllation between the two at the metro level.
// ---
// Alex Weinberg, November 2019
// ***************************************************

*ssc install estout
use "../output/metro_merged_j2j_bds", clear

preserve
	keep if year>=2001
	sort year quarter
	collapse (sum) j2jhire mainb (mean) estabs_entry estabs, by(year quarter)
	gen j2jhire_rate 		= j2jhire/mainb
	collapse (mean) j2jhire_rate (mean) estabs_entry estabs, by(year)
	gen entry_rate 			= estabs_entry/ estabs
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	egen mean_j2jhire_rate 	= mean(j2jhire_rate)
	egen mean_entry_rate 	= mean(entry_rate)
	gen dj2j 				= j2jhire_rate - mean_j2jhire_rate
	gen dentry 				= entry_rate - mean_entry_rate
	egen dj2j_sd 			= sd(dj2j)
	egen dentry_sd 			= sd(dentry)
	gen j2j_scaled 			= dj2j/dj2j_sd
	gen entry_scaled 		= dentry/dentry_sd
	tsset year
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	tsset, clear
	export delimited j2jhire_rate entry_rate using "../../../2_MATLAB_code/Input_data/j2j_entry_matlab.csv", replace
restore

keep if year == 2006 | year == 2009

// J2J
gen j2jrate_quarterly = j2jhire / mainb
bysort msa year: egen j2jrate = mean(j2jrate_quarterly) // mean of quarterly rate
drop if quarter > 1

sort msa year
by msa: gen j2jrate_T = j2jrate[_n+1] 
by msa: gen j2jhire_T = j2jhire[_n+1] 

gen pctchange_j2j = (j2jrate_T - j2jrate) / j2jrate
gen dlog_j2j = log(j2jhire_T) - log(j2jhire)

// New Firms

gen newfirm_rate = estabs_entry / estabs

sort msa year
by msa: gen newfirm_rate_T = newfirm_rate[_n+1] 
by msa: gen estabs_entry_T = estabs_entry[_n+1] 

gen pctchange_newfirm = (newfirm_rate_T - newfirm_rate) / newfirm_rate
gen dlog_newfirm = log(estabs_entry_T) - log(estabs_entry)

drop if year > 2006
sort msa state pctchange_j2j dlog_j2j pctchange_newfirm dlog_newfirm

// REGRESSION ---------------------------------------

// label var dlog_j2j "Delta Log J2J"
// label var dlog_newfirm "Delta New Firm"
// label var pctchange_j2j "Pct. Change J2J"
// label var pctchange_newfirm "Pct. Change New Firm"

// qui reg dlog_j2j dlog_newfirm 
// estimates store m1, title("DLog Baseline")
// estadd local StateFE "No" , replace

// qui reg dlog_j2j dlog_newfirm i.state
// estimates store m2, title("DLog State FE")
// estadd local StateFE "Yes" , replace

// qui reg pctchange_j2j pctchange_newfirm 
// estimates store m3, title("Pct. Change Baseline")
// estadd local StateFE "No" , replace

// qui reg pctchange_j2j pctchange_newfirm i.state
// estimates store m4, title("Pct. Change State FE")
// estadd local StateFE "Yes" , replace

// esttab m1 m2 m3 m4, label drop(*state) star se r2 s(StateFE N, label("State FE"))

// EXPORT DATA FOR FINAL PLOT --------------------------------------------------

bysort state: gen nval = _N

// Export Level Changes in Rates

gen rate_change_j2j      	= j2jrate_T - j2jrate
gen rate_change_estabs 		= newfirm_rate_T - newfirm_rate

export delimited rate_change_j2j rate_change_estabs nval using "../../../2_MATLAB_code/Input_data/matlab_ratechange_j2j_newfirm.csv", replace

reg rate_change_j2j rate_change_estabs
predict unwtd_level
reg rate_change_j2j rate_change_estabs [aweight=mainb]
predict wtd_level

export delimited unwtd_level wtd_level using "../../../2_MATLAB_code/Input_data/predicted_ratechange_j2j.csv", replace

