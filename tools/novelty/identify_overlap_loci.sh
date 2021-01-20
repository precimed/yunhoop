#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script identifies overlapping loci across multiple loci files (which
# may requires manual cleaning of certain larger groups of specific loci),
# but WITHOUT ANY WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 2 ]
then
  echo "Usage:      sh identify_overlap_loci.sh list_of_loci_files outfile"
  echo "Arguments:  list_of_loci_files - file that includes paths to loci files with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP|locusnum_LeadSNP]) and tags"
  echo "            outfile - out file that contains shared loci groups"
  echo "Example:    sh identify_overlap_loci.sh mood_psych_conj005_loci.txt mood_psych_conj005_shared_loci.txt"
  exit 0
fi

list_of_loci_files=$1
outfile=$2

files=`awk '{print $1}' $list_of_loci_files | paste -s -d " "`
rm -f $outfile
cat $files | awk '{print $1,$2,$3}' | sort | uniq | while read line; do
    chr=`echo $line | cut -d' ' -f1`
    minbp=`echo $line | cut -d' ' -f2`
    maxbp=`echo $line | cut -d' ' -f3`
    num=`awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($2>=minbp && $2<=maxbp || $3>=minbp && $3<=maxbp || $2<minbp && $3>maxbp)' $files | wc -l`
    if [ $num -gt 1 ]; then
        echo $chr $minbp $maxbp
        awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($2>=minbp && $2<=maxbp || $3>=minbp && $3<=maxbp || $2<minbp && $3>maxbp) {print FILENAME,$1,$2,$3,$4}' $files | sed 's/ $//' >> $outfile
        echo "" >> $outfile
    fi
done
cat $list_of_loci_files | while read line; do
    file=`echo $line | cut -d' ' -f1`
    tag=`echo $line | cut -d' ' -f2`
    sed -i "s|$file|$tag|" $outfile
done

#join lines of each group
sed -i ':a;N;/\n$/!s/\n/; /;ta;P;d' $outfile
#remove duplicates and split into groups again
cat $outfile | sort | uniq | awk '{ print length, $0 }' | sort -n -r | cut -d' ' -f2- | sed 's/$/\n/' | sed 's/; /\n/g' > $outfile.tmp
mv $outfile.tmp  $outfile

echo "Please manually remove larger groups consist of below lines from $(basename $outfile):"
echo "---------------------------------------------------"
cat $outfile | grep -v '^$' | sort | uniq -c | sort | awk '$1>1 {print $2,$3,$4,$5,$6}'
echo "---------------------------------------------------"
