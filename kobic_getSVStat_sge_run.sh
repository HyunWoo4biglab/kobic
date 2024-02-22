#!/bin/bash
#$ -V
#$ -cwd
#$ -q kbb.q
#$ -N hw_py
#$ -e /BiO/scratch/users/minhak/hyunwoo/dnsv_analysis/rare_disease/log/
#$ -o /BiO/scratch/users/minhak/hyunwoo/dnsv_analysis/rare_disease/log/
#$ -pe pe_slots 1
#$ -M barbaric0096@gmail.com
#$ -m eas

## load environment
##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/BiO/apps/ETCHING-1.4.2a/lib/

##excuting command
python /BiO/scratch/users/minhak/hyunwoo/src/kobic_getSVStat_240222.py -i /BiO/scratch/users/minhak/hyunwoo/data/meta_data/samples_info_fastq_total_new_symbolic_240222.map -o /BiO/scratch/users/minhak/hyunwoo/dnsv_analysis/rare_disease/ -r /BiO/scratch/users/minhak/hyunwoo/dnsv_analysis/rare_disease/kobic_getSVStat_currentlyCompleted.txt -t 30
