#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to compare latest loci contained in a new file against
# previous ones contained in an old file, and identify and ouput overlapping
# and novel loci WITHOUT WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley
# Contributors: Guy, Kevin, Shahram, Naz, Weiqiu, et al.

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 3 ]; then
  echo "Usage:     sh identify_overlap_and_novel_loci.sh newlocifile oldlocifile outprefix"
  echo "Arguments: newlocifile - file that contains new loci with 3 or 4 columns (CHR, MinBP, MaxBP, [LeadSNP])"
  echo "           oldlocifile - file that contains old loci with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP])"
  echo "           outprefix   - prefix of the output files that contain overlapping and novel loci"
  echo "Example:   sh identify_overlap_and_novel_loci.sh mood_bip_loci_cond_001.csv mood_gwas_loci.csv mood_bip_cond_001"
  exit 0
fi

#-------------------------------------------------------------------------#

if [ ! -f $1 ]; then
    echo ""
    echo "NOTE: $1 does't exist"
    exit 1
fi
if [ ! -f $2 ]; then
    echo ""
    echo "NOTE: $2 does't exist"
    exit 1
fi

if [ `grep ',' $1 | wc -l` -gt 0 ]; then
    sed 's/,/ /g' $1 > $1.tmp
    mv $1.tmp $1
fi
if [ `grep ',' $2 | wc -l` -gt 0 ]; then
    sed 's/,/ /g' $2 > $2.tmp
    mv $2.tmp $2
fi

if [ `head -n1 $1 | awk '{print NF}'` -ne 3 ] && [ `head -n1 $1 | awk '{print NF}'` -ne 4 ]; then
    echo ""
    echo "NOTE: $1 should contains 3 or 4 columns"
    exit 1
fi
if [ `head -n1 $2 | awk '{print NF}'` -ne 3 ] && [ `head -n1 $2 | awk '{print NF}'` -ne 4 ]; then
    echo ""
    echo "NOTE: $2 should contains 3 or 4 columns"
    exit 1
fi

tr -d '\r' < $1 > $1.tmp
mv $1.tmp $1

tr -d '\r' < $2 > $2.tmp
mv $2.tmp $2

files=$1" "$2
outprefix=$3
rm -f $outprefix.tmp
cat $files | awk '{print $1,$2,$3}' | sort | uniq | while read line; do
    chr=`echo $line | cut -d' ' -f1`
    minbp=`echo $line | cut -d' ' -f2`
    maxbp=`echo $line | cut -d' ' -f3`
    num=`awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($2>=minbp && $2<=maxbp || $3>=minbp && $3<=maxbp || $2<minbp && $3>maxbp)' $files | wc -l`
    if [ $num -gt 1 ]; then
        awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && ($2>=minbp && $2<=maxbp || $3>=minbp && $3<=maxbp || $2<minbp && $3>maxbp) {print FILENAME,$1,$2,$3,$4}' $files | sed 's/ $//' >> $outprefix.tmp
        echo "" >> $outprefix.tmp
    fi
done

sed ':a;N;/\n$/!s/\n/; /;ta;P;d' $outprefix.tmp > $outprefix.tmp2
mv $outprefix.tmp2 $outprefix.tmp

cat $outprefix.tmp | sort | uniq | awk '{ print length, $0 }' | sort -n -r | cut -d' ' -f2- | sed 's/$/\n/' | sed 's/; /\n/g' > ${outprefix}_group_loci.txt

grep $1 ${outprefix}_group_loci.txt | cut -d' ' -f2- | sed 's/ /	/g' | sort | uniq | sed 's/n$//g' > $outprefix.tmp
cat $1 | sed 's/ /	/g' | sort > $1.tmp
diff $1.tmp $outprefix.tmp | grep \< | cut -d' ' -f2- | sort -n > ${outprefix}_novel_loci.txt
sort -n $outprefix.tmp > ${outprefix}_overlap_loci.txt
rm -f $outprefix.tmp $1.tmp
rm -f ${outprefix}_group_loci.txt

echo "Done"
echo "Overlapping loci are saved in ${outprefix}_overlap_loci.txt"
echo "Novel loci are saved in ${outprefix}_novel_loci.txt"
