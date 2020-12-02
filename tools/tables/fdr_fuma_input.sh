#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates a snp list for FUMA analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 4 ]; then
  echo "Usage:     sh fdr_fuma_input.sh cfdr_clump_snp_file fdr r2 outfile"
  echo "Arguments: cfdr_clump_snp_file - file that contains cfdr clumping snp info"
  echo "           fdr - FDR filter for selecting snps"
  echo "           r2 - r2 filter for selecting snps"
  echo "           outfile - output file including selected snp list"
  echo "Example:   sh fdr_fuma_input.sh PGC_BIP_2016_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.snps.csv 0.1 0.6 BIP_vs_MOOD_input_snps.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

cfdr_clump_snp_file=$1
fdr=$2
r2=$3
outfile=$4

echo 'snp	pval' > $outfile
tail -n +2 $cfdr_clump_snp_file | awk -v fdr=$fdr -v r2=$r2 'NF==11 && $7>=r2 && $11<fdr {print $6,$11/1000000000}' OFS="\t" | sort | uniq >> $outfile
