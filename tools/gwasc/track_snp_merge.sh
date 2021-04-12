#!/bin/bash
#--------------------------- Description ---------------------------------#

# 
# 

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 3 ]; then
  echo "Usage: sh track_snp_merge.sh snp_merge_map snp_checklist outfile"
  echo "Arguments: snp_merge_map - map compiled from refsnp-merged.json"
  echo "           snp_checklist - file including the list of snps to check"
  echo "           outfile - output file"
  echo "Example: sh track_snp_merge.sh dbsnp_b153_merge_map.txt dbsnp_b153_merge_map_A_001.txt dbsnp_b153_merge_map_B_001.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

snp_merge_map=$1
snp_checklist=$2
outfile=$3

awk -F '\t' 'NR==FNR {D[$1]++; next} !($2 in D)' $snp_merge_map $snp_checklist | cut -f1-2 > $outfile
awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' $snp_merge_map $snp_checklist | cut -f1-2 > $snp_checklist.tmp
loop_flag='Y'
while [ "$loop_flag" = "Y" ]; do
    loop_flag='N'
    rm -f $snp_checklist.tmp2
    cat $snp_checklist.tmp | while read line; do
        merged_snp=`echo "$line" | cut -f1`
        merged_into_snp=`echo "$line" | cut -f2`
        grep "^${merged_into_snp}	" $snp_merge_map | cut -f2 | while read merged_into_snp2; do
            echo "${merged_snp}	$merged_into_snp2" >> $snp_checklist.tmp2
        done
    done
    if [ -s $snp_checklist.tmp2 ]; then
        awk -F '\t' 'NR==FNR {D[$1]++; next} !($2 in D)' $snp_merge_map $snp_checklist.tmp2 >> $outfile
        awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' $snp_merge_map $snp_checklist.tmp2 > $snp_checklist.tmp
        if [ -s $snp_checklist.tmp ]; then
            loop_flag='Y'
        fi
    fi
done
rm -f $snp_checklist.tmp*
