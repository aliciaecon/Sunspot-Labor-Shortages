#!/bin/sh

#====================#
# Get J2J Metro Data
#====================#

baseURL="https://lehd.ces.census.gov/data/j2j/R2019Q1/metro/j2j/"
manifest="https://lehd.ces.census.gov/data/j2j/R2019Q1/metro/j2j/j2j_metro_manifest.txt"

wget -O ../input/metromanifest.txt $manifest

echo "Downloading Metro-Level Datasets, will take a long time..."
while read ll; do
	datalink="$baseURL$ll"

	# download
	wget -P ../input $datalink
done < ../input/metromanifest.txt

# unzip .csv.gz files
gzip -dfr ../input

# Download geography labels
metro_label=https://lehd.ces.census.gov/data/schema/j2j_latest/label_geography_metro.csv
wget -O ../input/metro_label.csv $metro_label


#====================#
# Get BDS Data
#====================#

# no 2015-2016 data because they were updating the system that year
legacyURL=http://www2.census.gov/ces/bds/firm/bds_f_msa_release.csv
wget -O ../input/metroBDS.csv $legacyURL

#====================#
# Get MSA Codes
#====================#
# need this to create crosswalk
msa18=https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2018/delineation-files/list1_Sep_2018.xls
wget -O ../input/msa18.xls $msa18

msa9=https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2009/historical-delineation-files/list1.txt
wget -O ../input/msa9.txt $msa9
