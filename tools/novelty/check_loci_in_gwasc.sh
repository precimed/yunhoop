#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to check whether the loci reported by fdr analysis have
# been registered in gwascatalog with respect to specific phenotype, in which 
# case the loci are considered unnovel WITHOUT WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Weiqiu Cheng, Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 4 ]; then
  echo "Usage: sh check_loci_in_gwasc.sh fdr_clump_snp_file fuma_gwascatalog_file keywords"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           fuma_gwascatalog_file - fuma gwascatalog file"
  echo "           keywords - keywords with respect to specific phenotype, can be multiple words combined with & and multiple patterns delimited with | within quotes"
  echo "           outfile - output file including hits of gwascatalog"
  echo "Example: sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.txt 'high & density & lipoprotein | hdl & cholesterol' hdl_gwasc.txt"
  exit 0
fi

#-------------------------------------------------------------------------#

fdr_snps=$1
gwasc=$2
keywords=`echo $3 | sed "s/ //g"`
outfile=$4

sort -s -k6,6 $fdr_snps > $fdr_snps.sorted
sort -s -k5,5 $gwasc > $gwasc.sorted

join -1 6 -2 5 $fdr_snps.sorted $gwasc.sorted -t '	' | cut -f1-3,8,22-23 | sort -n -k2,2 | awk -F '\t' '{print $2,$3,$5,$6}' OFS='\t' > $fdr_snps.gwasc.sorted

n=1
for k in `echo $keywords | awk -F '|' '{for(i=1;i<=NF;i++){printf " %s", $i}}'`; do
    cp $fdr_snps.gwasc.sorted $fdr_snps.gwasc.$n.sorted
    for w in `echo $k | awk -F '&' '{for(i=1;i<=NF;i++){printf " %s", $i}}'`; do
        grep -i $w $fdr_snps.gwasc.$n.sorted > $fdr_snps.gwasc.$n.sorted.tmp
        mv $fdr_snps.gwasc.$n.sorted.tmp $fdr_snps.gwasc.$n.sorted
    done
    n=$((n+1))
done

cat $fdr_snps.gwasc.*.sorted | sort | uniq | sort -n -k1,1 > $outfile

rm -f ${fdr_snps}*.sorted $gwasc.sorted
