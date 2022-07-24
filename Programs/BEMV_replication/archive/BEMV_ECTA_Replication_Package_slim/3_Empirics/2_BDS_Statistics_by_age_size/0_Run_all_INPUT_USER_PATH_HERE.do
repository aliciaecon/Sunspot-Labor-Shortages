/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* THIS IS THE ONLY LINE USER NEEDS TO EDIT */
// Input path for the replication zip file:
gl root     	= "C:/Users/Simon Mongey/Dropbox/34_Bilal Engbom Mongey Violante/3_NewCode/Replication121221"
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* BELOW CODE RUNS ALL */
gl base 		= "$root/3_Empirics/2_BDS_Statistics_by_age_size"
gl working 		= "$root/3_Empirics/2_BDS_Statistics_by_age_size/working"
cd "$base"
do "$base/1_BEMV_Import_BDS_J2J"
do "$base/2_BEMV_Data_For_Paper"

