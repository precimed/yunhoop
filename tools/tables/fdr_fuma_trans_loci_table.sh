#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to identify transdiagnostic loci shared across multiple
# traits from pleioFDR/FUMA analysis WITHOUT WARRANTY.

# Authors: Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 2 ]
then
  echo "Usage:      sh fdr_fuma_trans_loci.sh list_of_loci_files outfile"
  echo "Arguments:  list_of_loci_files - file that includes paths to shared loci file, and fdr_fuma_snp_table, make sure that absolute or relative path to be specified"
  echo "            outfile - output file that contains transdiagnostic loci fdr/fuma results"
  echo "Example:    sh fdr_fuma_trans_loci.sh mood_psych_trans_conj005.txt mood_psych_conj005_trans_loci.txt"
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
    cut -f1-4,8-10 $trait1 | sort -s -k3,3 > $trait1.tmp
    cut -f1-4,8-10 $trait2 | sort -s -k3,3 > $trait2.tmp
    join -1 3 -2 3 $trait1.tmp $trait2.tmp | cut -d' ' -f1,3-5,11 | awk '{if($4>=$5) print $1,$2,$3,$4,$5,$4; else print $1,$2,$3,$4,$5,$5}' > snp_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$2==chr && $3>=minbp && $3<=maxbp' snp_results.csv | cut -d' ' -f1,6 | sort -n -k2,2 | head -n1
    rm -f $trait1.tmp $trait2.tmp $trait3.tmp snp_results.csv
}

function identify_low_max_snp_among_three_traits {
    trait1=$1
    trait2=$2
    trait3=$3
    chr=$4
    minbp=$5
    maxbp=$6
    cut -f1-4,8-10 $trait1 | sort -s -k3,3 > $trait1.tmp
    cut -f1-4,8-10 $trait2 | sort -s -k3,3 > $trait2.tmp
    cut -f1-4,8-10 $trait3 | sort -s -k3,3 > $trait3.tmp
    join -1 3 -2 3 $trait1.tmp $trait2.tmp > 2.tmp
    join -1 1 -2 3 2.tmp $trait3.tmp | cut -d' ' -f1,3-5,11,17 | awk '{if($4>=$5 && $4>=$6) print $1,$2,$3,$4,$5,$6,$4; else if($5>=$4 && $5>=$6) print $1,$2,$3,$4,$5,$6,$5; else print $1,$2,$3,$4,$5,$6,$6}' > snp_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$2==chr && $3>=minbp && $3<=maxbp' snp_results.csv | cut -d' ' -f1,7 | sort -n -k2,2 | head -n1
    rm -f $trait1.tmp $trait2.tmp $trait3.tmp 2.tmp snp_results.csv
}

function identify_low_max_snp_among_four_traits {
    trait1=$1
    trait2=$2
    trait3=$3
    trait4=$4
    chr=$5
    minbp=$6
    maxbp=$7
    cut -f1-4,8-10 $trait1 | sort -s -k3,3 > $trait1.tmp
    cut -f1-4,8-10 $trait2 | sort -s -k3,3 > $trait2.tmp
    cut -f1-4,8-10 $trait3 | sort -s -k3,3 > $trait3.tmp
    cut -f1-4,8-10 $trait4 | sort -s -k3,3 > $trait4.tmp 
    join -1 3 -2 3 $trait1.tmp $trait2.tmp > 2.tmp
    join -1 1 -2 3 2.tmp $trait3.tmp > 3.tmp
    join -1 1 -2 3 3.tmp $trait4.tmp | cut -d' ' -f1,3-5,11,17,23 | awk '{if($4>=$5 && $4>=$6 && $4>=$7) print $1,$2,$3,$4,$5,$6,$7,$4; else if($5>=$4 && $5>=$6 && $5>=$7) print $1,$2,$3,$4,$5,$6,$7,$5; else if($6>=$4 && $6>=$5 && $6>=$7) print $1,$2,$3,$4,$5,$6,$7,$6; else print $1,$2,$3,$4,$5,$6,$7,$7}' > snp_results.csv
    awk -v chr=$chr -v minbp=$minbp -v maxbp=$maxbp '$2==chr && $3>=minbp && $3<=maxbp' snp_results.csv | cut -d' ' -f1,8 | sort -n -k2,2 | head -n1
    rm -f $trait1.tmp $trait2.tmp $trait3.tmp $trait4.tmp 2.tmp 3.tmp snp_results.csv
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
    snpfuma1=snpfuma_$trait1
    snpfuma2=snpfuma_$trait2

    if [ $n_traits -ge 3 ]; then
        trait3=`echo $traits | cut -d' ' -f3`
        snpfuma3=snpfuma_$trait3
    fi

    if [ $n_traits -ge 4 ]; then
        trait4=`echo $traits | cut -d' ' -f4`
        snpfuma4=snpfuma_$trait4
    fi

    if [ $n_traits -eq 2 ]; then
        line=`identify_low_max_snp_among_two_traits ${!snpfuma1} ${!snpfuma2} $chr $min_overlap_bp $max_overlap_bp`
    elif [ $n_traits -eq 3 ]; then
        line=`identify_low_max_snp_among_three_traits ${!snpfuma1} ${!snpfuma2} ${!snpfuma3} $chr $min_overlap_bp $max_overlap_bp`
    elif [ $n_traits -eq 4 ]; then
        line=`identify_low_max_snp_among_four_traits ${!snpfuma1} ${!snpfuma2} ${!snpfuma3} ${!snpfuma4} $chr $min_overlap_bp $max_overlap_bp`
    fi
    if [ "$line" != "" ]; then
        rs=`echo $line | cut -d' ' -f1`
        lowmaxfdr=`echo $line | cut -d' ' -f2`
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
    else
        echo "no trans loci"
    fi
done
rm -f whatever-*.txt
