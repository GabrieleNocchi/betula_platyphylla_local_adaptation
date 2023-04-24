#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=24:0:0




module load python/3.8.5
source /data/home/btx388/RONA/bin/activate
python /data/home/btx388/RONA/bin/pyRona lfmm -pc current.txt -fc future_ssp3_2100.txt -out output.pdf -P 0.2 -assoc assoc_table.csv -geno birch.lfmm -covars 1
