_________________________________________________________________________________________
OVERVIEW:
_________________________________________________________________________________________
README file for 'Firm and Worker Dynamics in a Frictional Labor Market'
- Adrien Bilal, Niklas Engbom, Simon Mongey and Giovanni Violante
- Econometrica
- Written: January, 2022
_________________________________________________________________________________________
SUMMARY:
_________________________________________________________________________________________
--- This ZIP directory contains the following directories:
	--- ~/1_Draft
	--- ~/2_MATLAB_code
	--- ~/3_Empirics
--- We describe each of these below
_________________________________________________________________________________________
~/1_Draft
_________________________________________________________________________________________
--- Contains tex files for manuscript and online appendix.
--- Opening the tex files will allow a user to disambiguate the origin of any figure or table file
--- Compiling:
	--- ~1_Draft/BEMV_ECTA_Manuscript.tex
		--- Compiles jointly the manuscript and print appendix
	--- ~1_Draft/BEMV_ECTA_ONLINE_appendix.tex
		--- Compiles the online appendix
	--- The manuscript and appendix were compiled using PDFLaTeX in WinEdt v10.3
	--- Other tex files are inputs to these
	--- Both tex files can be compiled immediately after download
--- Also contains corresponding pdf files
--- Also contains figures that are drawn using tikz and don't need data or model
_________________________________________________________________________________________
~/2_MATLAB_code 
*** RUN AFTER 3_Empirics ***
_________________________________________________________________________________________
-1- The zip file comes with the following directories populated with 
	all figures and tables in the paper:
	--- Created_figure_files
	--- Created_table_files
	These are loaded into the tex document by the above tex files
	--- There is also a directory with outputs created by MATLAB
		during computations. 
		--- So that the replication package is not too big, it comes empty.
		--- It is called Created_mat_files
-2- To check replication of all results in the paper, first delete all files
	in the Created_* directories
-3- In MATLAB run ~/2_MATLAB_code/RUN_make_all.m
	--- This will populate all Created_* directories with 
		all figures and tables
	--- Then compile the tex files in ~/1_Draft which read in the objects from the 
		~/2_MATLAB_code/Created_* directories
	--- Doing so will replicate the printed manuscript and appendices.
-4- The Input_data directory contains:
	-1- Data produced in the 3_Empirics folder (see below)
		-A- From ~/3_Empriics/2_BDS_Statistics_by_age_size
			--- BDS_J2J_by_age.xls
			--- BDS_J2J_by_size.xls
		-B- From ~/3_Empriics/2_BDS_Statistics_by_age_size
			--- j2j_entry_matlab.xls	
			--- matlab_ratechange_j2j_newfirm.xls
			--- predicted_ratechange_j2j.xls
	-2- Output from Mongey_Violante_BLS_JOLTS_Microdata_Output.xls
		which is an appoved disclosure from BLS JOLTS/QCEW microdata
	-3- Targets_for_calibration.m, which is the list of target moments
		used in the production of estimation tables and figures
	-4- rdq028_Supplementary_Data which is supplementary data from Luttmer (ReStud 2011) 
		--- Download Supplementary data to Luttmer (2011)
		--- From: https://academic.oup.com/restud/article/78/3/1042/1566375?searchresult=1#supplementary-data
-5- The Input_mat_files directory contains: 
	--- Outcomes from the estimation of the model, which is not covered in the replication
	-1- Vector of estimated parameters
	-2- Jacobian of the objective function
_________________________________________________________________________________________
~/3_Empirics 
*** RUN BEFORE 2_MATLAB_code ***
_________________________________________________________________________________________
--- This directory contains all empirics

~/3_Empirics/1_Cross_MSA_Entry_EE
	--- Follow the enclosed README file.
	--- This generates Figure 13 in the paper
	
~/3_Empirics/2_BDS_Statistics_by_age_size
	--- Follow the enclosed README file
	--- This generates xls files that are then inputs to the MATLAB code
	--- These populate ~/2_MATLAB_code/Input_data
		--- BDS_J2J_by_age.xls
		--- BDS_J2J_by_size.xls
_________________________________________________________________________________________
	
	