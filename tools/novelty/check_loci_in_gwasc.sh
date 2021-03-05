#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to check whether the loci reported by fdr analysis have
# been registered in gwascatalog with respect to specific phenotype, in which 
# case the loci are considered to be not novel WITHOUT WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Weiqiu Cheng, Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 4 ]; then
  echo "Usage: sh check_loci_in_gwasc.sh fdr_clump_snp_file gwascatalog_file keywords outfile"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           gwascatalog_file - gwascatalog file"
  echo "           keywords - keywords with respect to specific phenotype, can be multiple words combined with & and multiple patterns delimited with | within quotes"
  echo "           outfile - output file including hits of gwascatalog"
  echo "Example: sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.tsv 'high & density & lipoprotein | hdl & cholesterol' hdl_gwasc.txt"
  exit 0
fi

#-------------------------------------------------------------------------#

fdr_clump_snp_file=$1
gwascatalog_file=$2
keywords=`echo $3 | sed "s/ //g"`
outfile=$4

cut -f1-2,6,8 $fdr_clump_snp_file | sort -s -k3,3 > $fdr_clump_snp_file.sorted
#FUMA GWAS Catalog
cut -f5,12-13 $gwascatalog_file | sort -s -k1,1 > $gwascatalog_file.sorted
#NHGRI-EBI GWAS Catalog
#cut -f7-8,22 $gwascatalog_file | awk -F '\t' '{print $3,$1,$2}' OFS='\t' | sort -s -k1,1 > $gwascatalog_file.sorted

join -1 3 -2 1 $fdr_clump_snp_file.sorted $gwascatalog_file.sorted -t '	' | sort -n -k2,2 | awk -F '\t' '{print $2,$3,$5,$6}' OFS='\t' > $fdr_clump_snp_file.gwasc.sorted

n=1
for k in `echo $keywords | awk -F '|' '{for(i=1;i<=NF;i++){printf " %s", $i}}'`; do
    cp $fdr_clump_snp_file.gwasc.sorted $fdr_clump_snp_file.gwasc.$n.sorted
    for w in `echo $k | awk -F '&' '{for(i=1;i<=NF;i++){printf " %s", $i}}'`; do
        grep -i $w $fdr_clump_snp_file.gwasc.$n.sorted > $fdr_clump_snp_file.gwasc.$n.sorted.tmp
        mv $fdr_clump_snp_file.gwasc.$n.sorted.tmp $fdr_clump_snp_file.gwasc.$n.sorted
    done
    n=$((n+1))
done

cat $fdr_clump_snp_file.gwasc.*.sorted | sort | uniq | sort -n -k1,1 > $outfile
rm -f ${fdr_clump_snp_file}*.sorted $gwascatalog_file.sorted
