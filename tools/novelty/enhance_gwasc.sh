#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to enhance NHGRI-EBI GWAS Catalog to make a clean file
# for novelty checking of loci WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -ne 2 ]; then
  echo "Usage: sh enhance_gwasc.sh gwasc outfile"
  echo "Arguments: gwasc - gwascatalog"
  echo "           outfile - output file"
  echo "Example: sh enhance_gwasc.sh gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv gwas_catalog_v1.0-associations_e100_r2021-02-25.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
outfile=$2

cut -f7-8,12-13,22 $gwasc | awk -F '\t' '{print $3,$4,$5,$1,$2}' OFS='\t' | tr -dc '\0-\177' | sed 's/&beta;/Beta /g' > $gwasc.txt

#filter irregular snps with empty chr
cut -f1-3 $gwasc.txt | awk -F '\t' '$1==""' | grep -v ":" | grep -i -v rs | grep -i -v chr | cut -f3 > $gwasc.tmp
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $gwasc.tmp $gwasc.txt > $gwasc.tmp2
mv $gwasc.tmp2 $gwasc.txt

#clean snps
cut -f1-3 $gwasc.txt | sed 's/\./:/' | sed 's/\s+//g' | sed 's/ x /;/gi' | sed 's/ ; /;/g' | sed 's/; /;/g' | sed 's/ ;/;/g' | sed 's/ , /,/g' | sed 's/, /,/g' | sed 's/ ,/,/g' | sed 's/,/,/g' | sed 's/ \/ /,/g' | sed 's/ \//,/g' | sed 's/\/ /,/g' | sed 's/\//,/g' > $gwasc.tmp1
cut -f4-5 $gwasc.txt > $gwasc.tmp2
for((i=0; i<=22; i++)); do
    sed "s/chr${i}: /chr${i}:/gi" $gwasc.tmp1 > $gwasc.tmp3
    sed "s/chr${i}_/chr${i}:/gi" $gwasc.tmp3 > $gwasc.tmp1
done
sed 's/chr://gi' $gwasc.tmp1 | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/del-//gi' | sed 's/		.*psy_/		/gi' > $gwasc.tmp3
paste -d'	' $gwasc.tmp3 $gwasc.tmp2 > $gwasc.tmp

#process lines with empty chr
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/ && $3!~/hg18/' | sed 's/	/|/g' > $gwasc.tmp1
rm -f $gwasc.tmp10

cut -d'|' -f3- $gwasc.tmp1 | while read line; do
    study=`echo $line | cut -d'|' -f2-3`
    echo $line | cut -d'|' -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | while read line; do
        snp=`echo $line | cut -d'_' -f1 | cut -d'-' -f1 | cut -d':' -f1-2 | sed 's/[a-zA-Z]$//'`
        if [ `echo $snp | grep ':' | wc -l` -gt 0 ]; then
            chr=`echo $snp | cut -d':' -f1`
            pos=`echo $snp | cut -d':' -f2`
            echo $chr'|'$pos'|'$snp'|'$study | sed 's/|/	/g' >> $gwasc.tmp10
        else
            echo '-|-|'$snp'|'$study | sed 's/|/	/g' >> $gwasc.tmp10
        fi
    done 
done
mv $gwasc.tmp10 $gwasc.tmp1

#process lines with single snp with chr
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/' $gwasc.tmp | tail -n +2 | cut -f1-5 > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f4-5 $gwasc.tmp2 > $gwasc.tmp23
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 > $gwasc.tmp2

#process lines with multiple snps
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/)' $gwasc.tmp | sed 's/	/|/g' > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo $line | cut -d'|' -f1 | sed 's/;/\n/g' > $gwasc.tmp31
    echo $line | cut -d'|' -f2 | sed 's/;/\n/g' > $gwasc.tmp32
    echo $line | cut -d'|' -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/[a-zA-Z]$//' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
    study=`echo $line | cut -d'|' -f4-5`
    n1=`cat $gwasc.tmp31 | wc -l`
    n3=`cat $gwasc.tmp33 | wc -l`
    if [ $n1 -eq 1 ] && [ `echo $line | cut -d'|' -f3 | grep ',' | wc -l` -gt 0 ]; then
        for((i=2; i<=$n3; i++)); do
            head -n1 $gwasc.tmp31 >> $gwasc.tmp31
            head -n1 $gwasc.tmp32 >> $gwasc.tmp32
        done
    fi
    if [ $n1 -lt $n3 ]; then
        sort -n -k1.3n $gwasc.tmp33 > $gwasc.tmp34
        mv $gwasc.tmp34 $gwasc.tmp33
    fi
    rm -f $gwasc.tmp34
    for ((i=1; i<=$n3; i++)); do
        echo $study | sed 's/|/	/g' >> $gwasc.tmp34
    done
    paste -d'	' $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp33 $gwasc.tmp34 | sed 's/^		/-	-	/' >> $gwasc.tmp30
    rm -f $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp33 $gwasc.tmp34
done
mv $gwasc.tmp30 $gwasc.tmp3

echo 'CHR	POS	SNP	STUDY	TRAIT' > $outfile
cat $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3 | sed 's/^X	/23	/' | sed 's/	X:/	23:/' | sed 's/^Y	/24	/' | sed 's/	Y:/	24:/' | sort | uniq | sort -n -k1,1 -k2,2 > $gwasc.tmp
awk -F '\t' '$1 ~ /^[0-9]+$/' $gwasc.tmp >> $outfile
awk -F '\t' '$1 !~ /^[0-9]+$/' $gwasc.tmp >> $outfile

rm -f $gwasc.tmp* $gwasc.txt
