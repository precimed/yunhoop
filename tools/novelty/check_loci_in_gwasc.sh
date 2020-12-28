#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script checks whether the loci reported by fdr analysis have been
# registered in gwascatalog with respect to specific phenotype, if so the
# loci are considered unnovel, but WITHOUT ANY WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Weiqiu Cheng, Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 3 ]; then
  echo "Usage: sh check_loci_in_gwasc.sh fdr_clump_snp_file fuma_gwascatalog_file keyword"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           fuma_gwascatalog_file - fuma gwascatalog file"
  echo "           keyword such as name or partial name of specified trait, can be multiple words within ''"
  echo "Example: sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.txt 'depress'"
  echo "         sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.txt 'major depressive disorder'"
  exit 0
fi

#-------------------------------------------------------------------------#

fdr_snps=$1
gwasc=$2
keyword=$3

sort -s -k6,6 $fdr_snps > $fdr_snps.sorted
sort -s -k5,5 $gwasc > $gwasc.sorted

join -1 6 -2 5 $fdr_snps.sorted $gwasc.sorted -t '	' | cut -f1-3,8,22 | sort -n -k2,2 | awk -F '\t' '{print $2,$3,$5}' OFS='\t' | grep -i "$keyword" | uniq

rm -f $fdr_snps.sorted $gwasc.sorted
