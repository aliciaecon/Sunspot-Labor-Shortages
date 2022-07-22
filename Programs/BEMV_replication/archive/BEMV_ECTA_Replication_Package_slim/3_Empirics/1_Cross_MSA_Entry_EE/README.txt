
FROM SRC: Run,...

metro_download.sh
--- Gets BDS and J2J data
--- ../input/metro_label.csv
--- ../input/metroBDS.csv
--- ../input/msa18.xls
--- ../input/msa9.txt

bds_clean.do
--- ../output/metroBDS

j2j_clean.do
--- ../output/j2j_metro

make_msa_crosswalk_2009_2018.do
--- ../output/crosswalk_09_18_v2

merge_metro.do
--- ../output/metro_merged_j2j_bds 

metro_compute_rates_plot.do
--- ../output/j2j_entry_matlab.csv
--- ../output/matlab_ratechange_j2j_newfirm.csv
--- ../output/predicted_ratechange_j2j.csv

plotTS_fmpv_bds_SM_AW.m
--- Figure_EE_Entry_data.eps
--- This figure appears in the manuscript.

