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
newlocifile=$1
oldlocifile=$2
outprefix=$3

if [ ! -f $newlocifile ]; then
    echo ""
    echo "NOTE: $newlocifile does't exist"
    exit 1
fi
if [ ! -f $oldlocifile ]; then
    echo ""
    echo "NOTE: $oldlocifile does't exist"
    exit 1
fi

if [ `grep ',' $newlocifile | wc -l` -gt 0 ]; then
    sed 's/,/ /g' $newlocifile > $newlocifile.tmp
    mv $newlocifile.tmp $newlocifile
fi
if [ `grep ',' $oldlocifile | wc -l` -gt 0 ]; then
    sed 's/,/ /g' $oldlocifile > $oldlocifile.tmp
    mv $oldlocifile.tmp $oldlocifile
fi

if [ `head -n1 $newlocifile | awk '{print NF}'` -ne 3 ] && [ `head -n1 $newlocifile | awk '{print NF}'` -ne 4 ]; then
    echo ""
    echo "NOTE: $newlocifile should contains 3 or 4 columns"
    exit 1
fi
if [ `head -n1 $oldlocifile | awk '{print NF}'` -ne 3 ] && [ `head -n1 $oldlocifile | awk '{print NF}'` -ne 4 ]; then
    echo ""
    echo "NOTE: $oldlocifile should contains 3 or 4 columns"
    exit 1
fi

tr -d '\r' < $newlocifile > $newlocifile.tmp
mv $newlocifile.tmp $newlocifile

tr -d '\r' < $oldlocifile > $oldlocifile.tmp
mv $oldlocifile.tmp $oldlocifile

files=$newlocifile" "$oldlocifile
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

grep $newlocifile ${outprefix}_group_loci.txt | cut -d' ' -f2- | sed 's/ /	/g' | sort | uniq | sed 's/n$//g' > $outprefix.tmp
cat $newlocifile | sed 's/ /	/g' | sort > $newlocifile.tmp
diff $newlocifile.tmp $outprefix.tmp | grep \< | cut -d' ' -f2- | sort -n > ${outprefix}_novel_loci.txt
sort -n $outprefix.tmp > ${outprefix}_overlap_loci.txt
rm -f $outprefix.tmp $newlocifile.tmp
rm -f ${outprefix}_group_loci.txt

echo "Overlapping loci are saved in ${outprefix}_overlap_loci.txt"
echo "Novel loci are saved in ${outprefix}_novel_loci.txt"
