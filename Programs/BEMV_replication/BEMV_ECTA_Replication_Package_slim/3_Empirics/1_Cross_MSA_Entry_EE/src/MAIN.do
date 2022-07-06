// ***************************************************
// This script runs all work on Metro data to produce
// input files for MATLAB code
// ***************************************************

// First run 1_metro_download in a bash / shell environment.
// This populates the 'input' folder (which is already populated in these replication materials)
	// Set directory of MAIN.do
cd "/Users/sebastianguarda/Dropbox/Bilal Engbom Mongey Violante/3_NewCode/Replication121221/3_Empirics/metro/src"

do 2_bds_clean
do 3_j2j_clean
do 4_make_msa_crosswalk_2009_2018
do 5_merge_metro
do 6_metro_compute_rates_plot
