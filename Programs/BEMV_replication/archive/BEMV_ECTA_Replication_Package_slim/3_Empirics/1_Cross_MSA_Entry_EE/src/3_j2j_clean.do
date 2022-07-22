// ***************************************************
// This script imports J2J csvs and merges them 
// 1. J2J from csv to dta
// 2. Append the J2Js from different years
// 3. Metro with the BDS metro data
// ---
// Alex Weinberg, November 2019
// ***************************************************


// CONVERT TO DTA ---------------------------------

// Metro Labels
import delimited using "../input/metro_label.csv", clear
// rename ïgeography geography
save "../input/metro_label", replace
clear all

// J2J files, takes a long time
local files : dir "../input" files "*.csv"

foreach file in `files' {
  // Only load J2J files
  if (substr("`file'",1,3) != "j2j") continue
  // drop .csv letters
  local outname = substr("`file'", 1, strlen("`file'") - 4)

  import delimited "../input/`file'", clear

  // keep only obs that include all race, sex, industry, etc.
  qui keep if agegrp=="A00"
  qui keep if education=="E0"
  qui keep if ethnicity=="A0"
  qui keep if firmage==0
  qui keep if firmsize==0
  qui keep if ownercode=="A00"
  qui keep if periodicity=="Q"
  qui keep if race=="A0"
  qui keep if sex==0
  qui keep if ind_level=="A"

  keep year quarter geography mainb j2jhire j2jsep
  drop if year > 2014 // bds metro data only until 2014
  save "../input/`outname'", replace
}


// APPEND ----------------------------------------

clear all
local files : dir "../input" files "*.dta"

foreach file in `files' {
  // only load j2j files
  if (substr("`file'",1,3) != "j2j") continue
  // vertically concatenate
  append using "../input/`file'"
}

// Merge Geography labels ---------------

merge m:1 geography using "../input/metro_label"
drop if _merge == 2 // just a handful of cities with no j2j data
drop _merge

// SAVE
rename geography msa18
label var mainb "Employment Beginning of Quarter"
save "../output/j2j_metro", replace


