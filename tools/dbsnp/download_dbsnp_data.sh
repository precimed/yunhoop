#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to 

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage: sh download_dbsnp_data.sh host out_folder"
  echo "Arguments: host - url to download build version related dbsnp data"
  echo "           out_folder - folder to place downloaded dbsnp data"
  echo "Example: sh download_dbsnp_data.sh https://ftp.ncbi.nih.gov/snp/redesign/archive/b153 $lf/ncbi/b153"
  exit 0
fi
#-------------------------------------------------------------------------#

host=$1
out_folder=$2

mkdir -p $out_folder/json
mkdir -p $out_folder/vcf
for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrMT merged withdrawn unsupported nosnppos other; do
    if [ ! -f $out_folder/json/refsnp-$i.json.bz2 ]; then
        echo $out_folder/json/refsnp-$i.json.bz2
        wget --continue --tries=0 $host/JSON/refsnp-$i.json.bz2 -O $out_folder/json/refsnp-$i.json.bz2
        wget --continue --tries=0 $host/JSON/refsnp-$i.json.bz2.md5 -O $out_folder/json/refsnp-$i.json.bz2.md5
        md5sum -c $out_folder/json/refsnp-$i.json.bz2.md5 > $out_folder/json/refsnp-checksum-$i.txt
        if [ "$i" = "merged" ] || [ "$i" = "nosnppos" ] || [ "$i" = "unsupported" ] || [ "$i" = "withdrawn" ] || [ "$i" = "other" ]; then
            bunzip2 -k $out_folder/json/refsnp-$i.json.bz2
        fi
    fi
done

for i in 25 38; do
    if [ ! -f $out_folder/vcf/GCF_000001405.${i}.gz ]; then
        wget --continue --tries=0 $host/VCF/GCF_000001405.${i}.gz -O $out_folder/vcf/GCF_000001405.${i}.gz
        wget --continue --tries=0 $host/VCF/GCF_000001405.${i}.gz.md5 -O $out_folder/vcf/GCF_000001405.${i}.gz.md5
        md5sum -c $out_folder/vcf/GCF_000001405.${i}.gz.md5 > $out_folder/vcf/gcf-checksum-$i.txt
    fi
done
