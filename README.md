# betula_platyphylla_local_adaptation
Code for:


Nocchi, G., Wang, J., Yang, L., Ding, J., Gao, Y., Buggs, R. J., & Wang, N. (2023). Genomic signals of local adaptation and hybridization in Asian white birch. Molecular Ecology, 32 (3), 595-612.

The script for the LFMM2 analysis is landgen1.R in the main directory.


The scripts for the RONA analysis are within the RONA directory:
- There are 7 directories within RONA, each corresponding to a different environmental variable.
- The input file in LFMM format containing SNP frequencies is birch.lfmm, zipped in birch.rar. This file is the same for all env. variables analyses.
- Each environemntal variable directory contains the environemntal data for that variable at present and in the future needed in pyRONA, file named accordingly.
- In each env. directory, there is an association table csv file. This file should have the p-values reported for the environemntal variable in question by LFMM, in the same order as found in the input birch.lfmm. The -P flag of pyRONA is used to tell pyRONA which SNP to include in the calculation of RONA, according to LFMM p-value. IMPORTANT, I have used an hacky solution here: to select significant SNPs I have converted p-values to q-values in the paper, which represent the FDR associated with a given p-value in our analysis. Therefore the association files that are found in these directories are custon made --> For each env. variable I gave SNPs that reported LFMM FDR < 0.1 for that variable a value of 0.1 in these assoc. table, while non-significant SNPs I gave them 0.9. In the pyRONA script (RONA.sh) I then use as threshold -P 0.2.  
- Outputs are in the files: RONA.sh.o[0-9]
- 1: AMT 2:MDR 3:ISO 4:MTWQ 5:AP 6:PS 7:PDQ


