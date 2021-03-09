#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to generate a snp list from pleioFDR output for FUMA
# analysis WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 4 ]; then
  echo "Usage:     sh fdr_fuma_snp_input.sh fdr_clump_snp_file fdr r2 outfile [exclude_region_list]"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           fdr - FDR filter for selecting snps"
  echo "           r2 - r2 filter for selecting snps"
  echo "           outfile - output file including selected snp list"
  echo "           exclude_region_list - regions with 'complex LD structure' to be excluded delimited by '#', eg:
                   MHC - 6 25119106 33854733; all traits. Given the size of the MHC it is probably sensible to exclude for all analyses
                   8p23 inversion - 8 7200000 12500000; psychiatric disorders and related traits
                   MAPT region - 17 43384864 44913631; neurological disorders e.g. Parkinson.s disease.
                   APOE region - 19 42000000 47000000; Alzheimer.s disease and autism spectrum disorder"
  echo "Example:   sh fdr_fuma_snp_input.sh conj.result.clump.snps.csv 0.1 0.6 TRAIT1_vs_TRAIT2_input_snps.txt '6 25119106 33854733 # 8 7200000 12500000'"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_clump_snp_file=$1
fdr=$2
r2=$3
outfile=$4

echo 'chr	bp	snp	pval' > $outfile
tail -n +2 $fdr_clump_snp_file | awk -v fdr=$fdr -v r2=$r2 'NF==11 && $7>=r2 && $11<fdr {print $2,$5,$6,$11/1000000000}' OFS="\t" | sort | uniq | sort -n -k1,1 -k2,2 >> $outfile

if [ $# -gt 4 ]; then
    exclude_region_list=$5
    rm -f $outfile.tmp
    echo $exclude_region_list | sed 's/#/\n/g' | while read i; do
        chr=`echo $i | awk '{print $1}'`
        minbp=`echo $i | awk '{print $2}'`
        maxbp=`echo $i | awk '{print $3}'`
        awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && $2>=minbp && $2<=maxbp' $outfile >> $outfile.tmp
    done
    sort $outfile.tmp | uniq | sort -n -k1,1 -k2,2 > $outfile.tmp2
    diff $outfile $outfile.tmp2 | grep \< | awk '{print $2,$3,$4,$5}' OFS='\t' > $outfile.tmp
    mv $outfile.tmp $outfile
fi
cut -f3-4 $outfile > $outfile.tmp2
mv $outfile.tmp2 $outfile
