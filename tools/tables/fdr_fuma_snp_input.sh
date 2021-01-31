#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates a snp list for FUMA analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 4 ]; then
  echo "Usage:     sh fdr_fuma_snp_input.sh fdr_clump_snp_file fdr r2 outfile"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           fdr - FDR filter for selecting snps"
  echo "           r2 - r2 filter for selecting snps"
  echo "           outfile - output file including selected snp list"
  echo "Example:   sh fdr_fuma_snp_input.sh conj.result.clump.snps.csv 0.1 0.6 trait1_vs_trait2_input_snps.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_clump_snp_file=$1
fdr=$2
r2=$3
outfile=$4

echo 'snp	pval' > $outfile
tail -n +2 $fdr_clump_snp_file | awk -v fdr=$fdr -v r2=$r2 'NF==11 && $7>=r2 && $11<fdr {print $6,$11/1000000000}' OFS="\t" | sort | uniq >> $outfile
