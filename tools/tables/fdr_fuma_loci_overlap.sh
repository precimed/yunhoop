#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to figure out overlaping loci of other traits for each
# specific trait across a group of FDR/FUMA analysis to fill table generated
# by script fdr_fuma_loci_table.sh TRAIT1_vs_TRAIT2_loci.txt WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh fdr_fuma_loci_overlap.sh fdr_fuma_loci_list outfile"
  echo "Arguments: fdr_fuma_loci_list - file that includes paths to TRAIT1_vs_TRAIT2_loci.txt files generated by script fdr_fuma_loci_table.sh"
  echo "           outfile - output file of overlaping loci"
  echo "Example:   sh fdr_fuma_loci_overlap.sh mood_psych_loci_conj005.txt mood_psych_conj005_shared_loci.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_fuma_loci_list=$1
outfile=$2

rm -f $(dirname $outfile)/trait_loci_list.txt
cat $fdr_fuma_loci_list | while read line; do
    fn=`echo $line | cut -d' ' -f1`
    trait=`echo $line | cut -d' ' -f2`
    head -n -1 $fn | awk '{print $2,$5,$6,$1"|"$3}' OFS='\t' > $(dirname $outfile)/tmp_trait_$trait.csv
    echo "$(dirname $outfile)/tmp_trait_$trait.csv	$trait" >> $(dirname $outfile)/trait_loci_list.txt
done

sh $(dirname $0)/../novelty/identify_overlap_loci.sh $(dirname $outfile)/trait_loci_list.txt $(dirname $outfile)/trait_shared_loci.txt
cut -d' ' -f1,5 $(dirname $outfile)/trait_shared_loci.txt | sed 's/|/ /g' | cut -d' ' -f1-2 | sed ':a;N;/\n$/!s/\n/; /;ta;P;d' > $outfile
echo ""
echo "Overlaping loci are included in $outfile"

rm -f $(dirname $outfile)/tmp_trait_*.csv $(dirname $outfile)/trait_loci_list.txt $(dirname $outfile)/trait_shared_loci.txt

cat $fdr_fuma_loci_list | while read line; do
    fn=`echo $line | cut -d' ' -f1`
    trait=`echo $line | cut -d' ' -f2`
    grep "${trait}" $outfile > ${outfile%.*}_$trait.tmp1
    grep "${trait}" $outfile | sed 's/; /\n/' | grep $trait > ${outfile%.*}_$trait.tmp2
    rm -f ${outfile%.*}_$trait.tmp
    paste -d'|' ${outfile%.*}_$trait.tmp2 ${outfile%.*}_$trait.tmp1 | while read i; do
        left=`echo $i | cut -d'|' -f1`
        right=`echo $i | cut -d'|' -f2 | sed "s/$left//" | sed 's/^; //' | sed 's/; ;/;/g' | sed 's/; $//'`
        echo $left $right | cut -d' ' -f2- >> ${outfile%.*}_$trait.tmp
    done
    sort -s -k1,1 ${outfile%.*}_$trait.tmp > ${outfile%.*}_$trait.tmp2
    head -n -1 $fn | tail -n +2 | cut -f1 | sort -s -k1,1 > $fn.tmp
    rm -f ${outfile%.*}_$trait.tmp
    join -1 1 -2 1 -a 1 $fn.tmp ${outfile%.*}_$trait.tmp2 | sort -n -k1,1 | awk '{a[$1]= a[$1]"; "$0} END {for (item in a ) print item, a[item]}' | sed 's/ ;//' | cut -d' ' -f2- | sort -n -k1,1 | while read i; do
        left=`echo $i | cut -d' ' -f1`
        right=`echo $i | cut -d' ' -f2- | sed "s/; $left/;/g"`
        if [ "$left" = "$right" ]; then
            echo $left >> ${outfile%.*}_$trait.tmp
        else
            echo $left $right >> ${outfile%.*}_$trait.tmp
        fi
    done
    
    echo "Overlaping_loci" > ${fn%.*}_$trait.txt
    cat ${outfile%.*}_$trait.tmp | cut -d' ' -f2- | awk '{if(NF==1) print ""; else print $0}' >> ${fn%.*}_$trait.txt
    paste -d'	' $fn ${fn%.*}_$trait.txt > ${fn%.*}_$trait.tmp
    mv ${fn%.*}_$trait.tmp ${fn%.*}.txt
    rm -f $fn.tmp ${outfile%.*}_$trait.tmp* ${fn%.*}_$trait.txt
    echo "See ${fn%.*}.txt for updated loci table"
done
