#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script identifies transdiagnostic loci shared across multiple traits
# from FDR/FUMA analysis, but WITHOUT ANY WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 2 ]
then
  echo "Usage:      sh fdr_fuma_trans_loci.sh list_of_loci_files outfile"
  echo "Arguments:  list_of_loci_files - file that includes paths to shared loci file, result.mat.csv of fdr, and fdr_fuma_snp_table, make sure that absolute or relative path to be specified"
  echo "            outfile - out file that contains transdiagnostic loci fdr/fuma results"
  echo "Example:    sh fdr_fuma_trans_loci.sh mood_psych_shared_loci_conj005.txt mood_psych_conj005_trans_loci.txt"
  exit 0
fi

#-------------------------------------------------------------------------#

list_of_loci_files=$1
outfile=$2
source $list_of_loci_files

function identify_low_max_snp_among_two_traits {
    trait1=$1
    trait2=$2
    chr=$3
    minbp=$4
    maxbp=$5
    paste -d'	' ${trait1} ${trait2} | cut -f1-6,12 | awk 'NF==7' | awk '{if($6>=$7) print $1,$2,$3,$4,$5,$6,$7,$6; else print $1,$2,$3,$4,$5,$6,$7,$7}' | awk '$8<1' > fdr_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && $3>=minbp && $3<=maxbp' fdr_results.csv | cut -d' ' -f1-5,8 | sort -n -k6,6 | head -n1
    rm -f fdr_results.csv
}

function identify_low_max_snp_among_three_traits {
    trait1=$1
    trait2=$2
    trait3=$3
    chr=$4
    minbp=$5
    maxbp=$6
    paste -d'	' $trait1 $trait2 $trait3 | cut -f1-6,12,18 | awk 'NF==8' | awk '{if($6>=$7 && $6>=$8) print $1,$2,$3,$4,$5,$6,$7,$8,$6; else if($7>=$6 && $7>=$8) print $1,$2,$3,$4,$5,$6,$7,$8,$7; else if($8>=$6 && $8>=$7) print $1,$2,$3,$4,$5,$6,$7,$8,$8}' | awk '$9<1' > fdr_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && $3>=minbp && $3<=maxbp' fdr_results.csv | cut -d' ' -f1-5,9 | sort -n -k6,6 | head -n1
    rm -f fdr_results.csv
}

function identify_low_max_snp_among_four_traits {
    trait1=$1
    trait2=$2
    trait3=$3
    trait4=$4
    chr=$5
    minbp=$6
    maxbp=$7
    paste -d'	' $trait1 $trait2 $trait3 $trait4 | cut -f1-6,12,18,24 | awk 'NF==9' | awk '{if($6>=$7 && $6>=$8 && $6>=$9) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$6; else if($7>=$6 && $7>=$8 && $7>=$9) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$7; else if($8>=$6 && $8>=$7 && $8>=$9) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$9}' | awk '$10<1' > fdr_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$1==chr && $3>=minbp && $3<=maxbp' fdr_results.csv | cut -d' ' -f1-5,10 | sort -n -k6,6 | head -n1
    rm -f fdr_results.csv
}

echo "Trait1	locusnum	CHR	bp_range	overlapping_bp_range	lowest_max_conjFDR	trans_Lead_SNP	trans_Lead_BP	conjFDR	A1	A2	PVAL_Trait1	PVAL_Trait2	Z_Trait1	Z_Trait2	BETA_Trait1	BETA_Trait2	cross-traitConcordant	Novel_in_Trait1	Novel_in_Trait2" > $outfile
awk -v RS= '{print > ("whatever-" NR ".txt")}' $shared_loci_file
n=`ls whatever-*.txt | wc -l`
for ((i=1;i<=$n;i++)); do
    traits=`awk '{print $1}' whatever-$i.txt | tr '\n' ' '`
    echo $traits
    n_traits=`cat whatever-$i.txt | wc -l`
    if [ $n_traits -lt 2 ]; then
        continue
    fi
    chr=`head -n1 whatever-$i.txt | awk '{print $2}'`
    min_overlap_bp=`sort -n -k3,3 whatever-$i.txt | tail -n1 | awk '{print $3}'`
    max_overlap_bp=`sort -n -k4,4 whatever-$i.txt | head -n1 | awk '{print $4}'`

    trait1=`echo $traits | cut -d' ' -f1`
    trait2=`echo $traits | cut -d' ' -f2`
    resfdr1=resfdr_$trait1
    resfdr2=resfdr_$trait2

    if [ $n_traits -ge 3 ]; then
        trait3=`echo $traits | cut -d' ' -f3`
        resfdr3=resfdr_$trait3
    fi

    if [ $n_traits -ge 4 ]; then
        trait4=`echo $traits | cut -d' ' -f4`
        resfdr4=resfdr_$trait4
    fi

    if [ $n_traits -eq 2 ]; then
        line=`identify_low_max_snp_among_two_traits ${!resfdr1} ${!resfdr2} $chr $min_overlap_bp $max_overlap_bp`
    elif [ $n_traits -eq 3 ]; then
        line=`identify_low_max_snp_among_three_traits ${!resfdr1} ${!resfdr2} ${!resfdr3} $chr $min_overlap_bp $max_overlap_bp`
    elif [ $n_traits -eq 4 ]; then
        line=`identify_low_max_snp_among_four_traits ${!resfdr1} ${!resfdr2} ${!resfdr3} ${!resfdr4} $chr $min_overlap_bp $max_overlap_bp`
    fi
    rs=`echo $line | cut -d' ' -f2`
    lowmaxfdr=`echo $line | cut -d' ' -f6`
    rm -f snpinfo.csv
    awk -v min_overlap_bp=$min_overlap_bp -v max_overlap_bp=$max_overlap_bp '{print $1,$2,$3,$4,min_overlap_bp,max_overlap_bp}' whatever-$i.txt | while read line; do
        trait=`echo $line | cut -d' ' -f1`
        minbp=`echo $line | cut -d' ' -f3`
        maxbp=`echo $line | cut -d' ' -f4`
        snpfuma=snpfuma_$trait
        locifuma=${!snpfuma/snps.txt/loci.txt}
        if [ ! -f $locifuma ]; then
            echo "$locifuma does not exist, please rename your fuma loci file to it"
            exit 1
        fi
        locusnum=`awk -v rs=$rs '$3==rs {print $1}' ${!snpfuma}`
        novelty1=`awk -v locusnum=$locusnum '$1==locusnum {print $26}' $locifuma`
        novelty2=`awk -v locusnum=$locusnum '$1==locusnum {print $27}' $locifuma`
        awk -v rs=$rs -v trait=$trait -v chr=$chr -v bp=$minbp-$maxbp -v overlapbp=$min_overlap_bp-$max_overlap_bp -v lowmaxfdr=$lowmaxfdr -v novelty1=$novelty1 -v novelty2=$novelty2 '$3==rs {print trait,$1,chr,bp,overlapbp,lowmaxfdr,rs,$4,$8,$9,$10,$24,$28,$25,$29,$26,$30,novelty1,novelty2}' ${!snpfuma} >> snpinfo.csv
    done
    n_trait1_pos=`awk '$16 >= 0' snpinfo.csv | wc -l`
    n_trait2_pos=`awk '$17 >= 0' snpinfo.csv | wc -l`
    n_trait1_neg=`awk '$16 < 0' snpinfo.csv | wc -l`
    n_trait2_neg=`awk '$17 < 0' snpinfo.csv | wc -l`
    if ([ $n_trait1_pos -eq 0 ] && [ $n_trait2_pos -eq 0 ]) || ([ $n_trait1_neg -eq 0 ] && [ $n_trait2_neg -eq 0 ]); then
        awk '{out=""; for(i=1;i<=17;i++){out=out" "$i}; print out,"TRUE",$18,$19}' snpinfo.csv | sed 's/^ //' | sed 's/ /	/g' >> $outfile
    else
        awk '{out=""; for(i=1;i<=17;i++){out=out" "$i}; print out,"FALSE",$18,$19}' snpinfo.csv | sed 's/^ //' | sed 's/ /	/g' >> $outfile
    fi
    rm -f snpinfo.csv
    echo "" >> $outfile
done
rm -f whatever-*.txt
