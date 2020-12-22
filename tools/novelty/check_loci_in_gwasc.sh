#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script checks whether the loci reported by cfdr analysis have been
# registered in gwascatalog with respect to specific phenotype, if so the
# loci are considered unnovel, but WITHOUT ANY WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Weiqiu Cheng, Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 3 ]; then
  echo "Usage: sh check_loci_in_gwasc.sh cfdr_clump_snps_file fuma_gwascatalog_file keyword"
  echo "Arguments: cfdr_clump_snps_file - file that contains snp info from clumping"
  echo "           fuma_gwascatalog_file - file that contains gwascatalog info from FUMA"
  echo "           keyword - keyword such as name or part name of investigated trait"
  echo "Example: sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.txt depress"
  exit 0
fi

#-------------------------------------------------------------------------#

cfdr_snps=$1
gwasc=$2
keyword=$3

sort -s -k6,6 $cfdr_snps > $cfdr_snps.sorted
sort -s -k5,5 $gwasc > $gwasc.sorted

join -1 6 -2 5 $cfdr_snps.sorted $gwasc.sorted -t '	' | cut -f1-3,8,22 | sort -n -k2,2 | awk -F '\t' '{print $2,$3,$5}' OFS='\t' | grep -i $keyword | uniq

rm -f $cfdr_snps.sorted $gwasc.sorted
