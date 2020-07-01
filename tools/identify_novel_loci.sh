#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script compares latest loci contained in a new file against previous
# ones contained in an old file, and identifies and ouputs the novel loci.

# Author: Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO
#-------------------------------------------------------------------------#

if [ $# -lt 3 ]; then
  echo "Usage:     sh identify_novel_loci.sh newlocifile oldlocifile outputfile"
  echo "Arguments: newlocifile - file that contains new loci with 4 columns (CHR, LEAD_SNP, MinBP, MaxBP)"
  echo "           oldlocifile - file that contains old loci with 4 columns (CHR, LEAD_SNP, MinBP, MaxBP)"
  echo "           outputfile  - output file that contains novel loci"
  echo "Example:   sh identify_novel_loci.sh mood_bip_loci_cond_001.csv mood_gwas_loci.csv novel_loci_cond_001_mood_bip.txt"
  exit 0
fi

files=$1" "$2
outputfile=$3
rm -f $outputfile
cat $files | cut -d$'\t' -f1,3-4 | sort | uniq | while read line; do
    chr=`echo $line | cut -d' ' -f1`
    minbp=`echo $line | cut -d' ' -f2`
    maxbp=`echo $line | cut -d' ' -f3`
    num=`awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($3>=minbp && $3<=maxbp || $4>=minbp && $4<=maxbp || $3<minbp && $4>maxbp) {print FILENAME,$1,$3,$4,$2}' $files | wc -l`
    if [ $num -gt 1 ]; then
        awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($3>=minbp && $3<=maxbp || $4>=minbp && $4<=maxbp || $3<minbp && $4>maxbp) {print FILENAME,$1,$3,$4,$2}' $files >> $outputfile.tmp
        echo "" >> $outputfile.tmp
    fi
done

sed -i ':a;N;/\n$/!s/\n/; /;ta;P;d' $outputfile.tmp
cat $outputfile.tmp | sort | uniq | awk '{ print length, $0 }' | sort -n -r | cut -d' ' -f2- | sed 's/$/\n/' | sed 's/; /\n/g' > $outputfile.tmp2
grep $1 $outputfile.tmp2 | awk '{print $2,$5,$3,$4}' OFS='\t' | sort | uniq  > $outputfile.tmp
sort $1 > $1.tmp
diff $1.tmp $outputfile.tmp | grep \< | cut -d' ' -f2- | sort -n > $outputfile
rm -f $outputfile.tmp* $1.tmp
