# betula_platyphylla_local_adaptation
Code for:


Nocchi, G., Wang, J., Yang, L., Ding, J., Gao, Y., Buggs, R. J., & Wang, N. (2023). Genomic signals of local adaptation and hybridization in Asian white birch. Molecular Ecology, 32 (3), 595-612.

The script for the LFMM2 analysis is landgen1.R in the main directory. 

The repository also include scripts to download environmental variables from WorldClim and assess their correlation, as well as scripts to analyze the structure of the input data set before GEA analysis.


Input files provided in the directory data. 
Elevation raster and original map/ped/nosex files not provided due to size restriction. 
Input provided in LFMM format for 71 Chinese nonhybrid B. platyphylla individuals. Geographic coordinates available in the file environmental_data.txt


The scripts for the RONA analysis are within the RONA directory:
- There are 7 directories within RONA, each corresponding to a different environmental variable.
  
- Each environmental variable directory contains the environemntal data for that variable at present and in the future which are needed in pyRONA; files are named accordingly.
  
- In each env. directory, there is an association table csv file. This file should have the SNPs p-values  for the environemntal variable in question reported by LFMM, given in the same SNP order as found in the input birch.lfmm. The -P flag of pyRONA is used to tell pyRONA which SNPs to include in the calculation of RONA, which are usually selected according to LFMM p-value.
  
- Outputs of each RONA run are in the files named: RONA.sh.o[0-9]* and output.pdf (which are plots created automatically by pyRONA).
  
- 1: AMT 2:MDR 3:ISO 4:MTWQ 5:AP 6:PS 7:PDQ


