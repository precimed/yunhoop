#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary loci table (cFDR+FUMA+SUMSTAT)
# from cFDR/FUMA analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -ne 3 ] && [ $# -ne 5 ] && [ $# -ne 9 ]; then
  echo "Usage:     sh fdr_fuma_loci_table.sh fdr_clump_loci_file fdr_fuma_snp_table outfile [trait1locifile trait2locifile] [fdr_snp_table fuma_gwasc keyword1 keyword2]"
  echo "Arguments: fdr_clump_loci_file - fdr clumping loci file"
  echo "           fdr_fuma_snp_table - file generated by script fdr_fuma_snp_table.sh"
  echo "           outfile - output file"
  echo "           trait1locifile - old loci file with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP]) of trait1, use '-' for none"
  echo "           trait2locifile - old loci file with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP]) of trait2, use '-' for none"
  echo "           fdr_snp_table - fdr clumping snp file, use '-' for none"
  echo "           fuma_gwasc - fuma gwas catalog file, use '-' for none"
  echo "           keyword1 - keyword of trait1, use '-' for none"
  echo "           keyword2 - keyword of trait2, use '-' for none"
  echo "Example:   sh fdr_fuma_loci_table.sh conj.result.clump.loci.csv trait1_vs_trait2_snps_conj.txt trait1_vs_trait2_loci.txt trait1_old_loci.csv trait2_old_loci.csv conj.result.clump.snps.csv gwascatalog.txt keyword1 keyword2"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_clump_loci_file=$1
fdr_fuma_snp_table=$2
outfile=$3

tag1=`basename $fdr_fuma_snp_table | cut -d'_' -f1`
tag2=`basename $fdr_fuma_snp_table | cut -d'_' -f3`
echo "locusnum	CHR	LEAD_SNP	LEAD_BP	MinBP	MaxBP	cFDR	non_effect_allele	effect_allele	nearestGene	dist	func	CADD	RDB	minChrState	commonChrState	${tag1}_PVAL	${tag1}_Z	${tag1}_BETA	${tag2}_PVAL	${tag2}_Z	${tag2}_BETA	Concordant_Effect" > $outfile

sort -s -k3,3 $fdr_clump_loci_file > $fdr_clump_loci_file.sorted
sort -s -k1,1 $fdr_fuma_snp_table > $fdr_fuma_snp_table.sorted

join -1 3 -2 1 -a 1 $fdr_clump_loci_file.sorted $fdr_fuma_snp_table.sorted | grep -v FDR | awk '{print $2,$3,$1,$4,$5,$6,$7,$12,$13,$17,$18,$19,$20,$21,$22,$23,$27,$28,$29,$30,$31,$32,$33}' OFS='\t' | sort -n -k2,2 -k4,4 >> $outfile

if [ $# -ge 5 ]; then
    trait1locifile=$4
    trait2locifile=$5
    tail -n +2 $outfile | awk '{print $2,$5,$6}' OFS='\t' > ${outfile%.*}.tmp

    if [ "$trait1locifile" != "-" ]; then
        awk '{print $1,$2,$3}' OFS='\t' $trait1locifile > ${trait1locifile%.*}.tmp
        sh $(dirname $0)/../novelty/identify_overlap_and_novel_loci.sh ${outfile%.*}.tmp ${trait1locifile%.*}.tmp ${outfile%.*}_${tag1}
        awk '{print $1"-"$2"-"$3}' ${outfile%.*}.tmp > ${outfile%.*}.tmp2
        awk '{print $1"-"$2"-"$3}' ${outfile%.*}_${tag1}_novel_loci.txt > ${outfile%.*}_${tag1}_novel_loci.tmp
        rm -f ${outfile%.*}_${tag1}_novel_loci.tmp2
        for i in `cat ${outfile%.*}.tmp2`; do
            if [ `grep "^$i$" ${outfile%.*}_${tag1}_novel_loci.tmp | wc -l` -gt 0 ]; then
                echo $i"	Yes" >> ${outfile%.*}_${tag1}_novel_loci.tmp2
            else
                echo $i"	No" >> ${outfile%.*}_${tag1}_novel_loci.tmp2
            fi
        done
    fi
    if [ "$trait2locifile" != "-" ]; then
        awk '{print $1,$2,$3}' OFS='\t' $trait2locifile > ${trait2locifile%.*}.tmp
        sh $(dirname $0)/../novelty/identify_overlap_and_novel_loci.sh ${outfile%.*}.tmp ${trait2locifile%.*}.tmp ${outfile%.*}_${tag2}
        awk '{print $1"-"$2"-"$3}' ${outfile%.*}.tmp > ${outfile%.*}.tmp2
        awk '{print $1"-"$2"-"$3}' ${outfile%.*}_${tag2}_novel_loci.txt > ${outfile%.*}_${tag2}_novel_loci.tmp
        rm -f ${outfile%.*}_${tag2}_novel_loci.tmp2
        for i in `cat ${outfile%.*}.tmp2`; do
            if [ `grep "^$i$" ${outfile%.*}_${tag2}_novel_loci.tmp | wc -l` -gt 0 ]; then
                echo $i"	Yes" >> ${outfile%.*}_${tag2}_novel_loci.tmp2
            else
                echo $i"	No" >> ${outfile%.*}_${tag2}_novel_loci.tmp2
            fi
        done
    fi
fi

if [ $# -eq 9 ]; then
    fdr_snp_table=$6
    fuma_gwasc=$7
    tail -n +2 $outfile | cut -f1 > ${outfile%.*}.tmp3
    if [ "$fdr_snp_table" != "-" ] && [ "$fuma_gwasc" != "-" ]; then
        keyword1=$8
        if [ "$keyword1" != "-" ]; then
            rm -f ${outfile%.*}_${tag1}_gwas.txt
            sh $(dirname $0)/../novelty/check_loci_in_gwasc.sh $fdr_snp_table $fuma_gwasc "$keyword1" > ${outfile%.*}_${tag1}_gwas.txt
            sh $(dirname $0)/../novelty/check_loci_in_gwasc.sh $fdr_snp_table $fuma_gwasc "$keyword1" | cut -f1 | uniq > ${outfile%.*}_${tag1}_novel_loci.tmp3
            rm -f ${outfile%.*}_${tag1}_novel_loci.tmp4
            for i in `cat ${outfile%.*}.tmp3`; do
                if [ `grep "^$i$" ${outfile%.*}_${tag1}_novel_loci.tmp3 | wc -l` -eq 0 ]; then
                    echo $i"	Yes" >> ${outfile%.*}_${tag1}_novel_loci.tmp4
                else
                    echo $i"	No" >> ${outfile%.*}_${tag1}_novel_loci.tmp4
                fi
            done
        fi
        keyword2=$9
        if [ "$keyword2" != "-" ]; then
            rm -f ${outfile%.*}_${tag2}_gwas.txt
            sh $(dirname $0)/../novelty/check_loci_in_gwasc.sh $fdr_snp_table $fuma_gwasc "$keyword2" > ${outfile%.*}_${tag2}_gwas.txt
            sh $(dirname $0)/../novelty/check_loci_in_gwasc.sh $fdr_snp_table $fuma_gwasc "$keyword2" | cut -f1 | uniq > ${outfile%.*}_${tag2}_novel_loci.tmp3
            rm -f ${outfile%.*}_${tag2}_novel_loci.tmp4
            for i in `cat ${outfile%.*}.tmp3`; do
                if [ `grep "^$i$" ${outfile%.*}_${tag2}_novel_loci.tmp3 | wc -l` -eq 0 ]; then
                    echo $i"	Yes" >> ${outfile%.*}_${tag2}_novel_loci.tmp4
                else
                    echo $i"	No" >> ${outfile%.*}_${tag2}_novel_loci.tmp4
                fi
            done
        fi
    fi
fi
if [ $# -ge 5 ]; then
    echo "locusnum	novel_in_novdb	novel_in_gwas" > ${outfile%.*}_${tag1}_novelty.txt
    if [ -f ${outfile%.*}_${tag1}_novel_loci.tmp2 ] && [ -f ${outfile%.*}_${tag1}_novel_loci.tmp4 ]; then
        paste ${outfile%.*}_${tag1}_novel_loci.tmp2 ${outfile%.*}_${tag1}_novel_loci.tmp4 | awk '{print $3,$2,$4}' OFS='\t' >> ${outfile%.*}_${tag1}_novelty.txt
    elif [ ! -f ${outfile%.*}_${tag1}_novel_loci.tmp2 ] && [ -f ${outfile%.*}_${tag1}_novel_loci.tmp4 ]; then
        awk '{print $1,"Yes",$2}' OFS='\t' ${outfile%.*}_${tag1}_novel_loci.tmp4 >> ${outfile%.*}_${tag1}_novelty.txt
    elif [ -f ${outfile%.*}_${tag1}_novel_loci.tmp2 ] && [ ! -f ${outfile%.*}_${tag1}_novel_loci.tmp4 ]; then
        awk '{print $0,"Yes"}' OFS='\t' ${outfile%.*}_${tag1}_novel_loci.tmp2 >> ${outfile%.*}_${tag1}_novelty.txt
    fi

    echo "locusnum	novel_in_novdb	novel_in_gwas" > ${outfile%.*}_${tag2}_novelty.txt
    if [ -f ${outfile%.*}_${tag2}_novel_loci.tmp2 ] && [ -f ${outfile%.*}_${tag2}_novel_loci.tmp4 ]; then
        paste ${outfile%.*}_${tag2}_novel_loci.tmp2 ${outfile%.*}_${tag2}_novel_loci.tmp4 | awk '{print $3,$2,$4}' OFS='\t' >> ${outfile%.*}_${tag2}_novelty.txt
    elif [ ! -f ${outfile%.*}_${tag2}_novel_loci.tmp2 ] && [ -f ${outfile%.*}_${tag2}_novel_loci.tmp4 ]; then
        awk '{print $1,"Yes",$2}' OFS='\t' ${outfile%.*}_${tag2}_novel_loci.tmp4 >> ${outfile%.*}_${tag2}_novelty.txt
    elif [ -f ${outfile%.*}_${tag2}_novel_loci.tmp2 ] && [ ! -f ${outfile%.*}_${tag2}_novel_loci.tmp4 ]; then
        awk '{print $0,"Yes"}' OFS='\t' ${outfile%.*}_${tag2}_novel_loci.tmp2 >> ${outfile%.*}_${tag2}_novelty.txt
    fi
    paste $outfile ${outfile%.*}_${tag1}_novelty.txt ${outfile%.*}_${tag2}_novelty.txt | cut -f1-23,25-26,28-29 | awk '{if($24=="Yes" && $25=="Yes") print $0,"Yes"; else print $0,"No"}' OFS='\t' | awk '{if($26=="Yes" && $27=="Yes") print $0,"Yes"; else print $0,"No"}' OFS='\t' | cut -f1-23,28-29 | sed "1s/No	No/Novel_in_${tag1}	Novel_in_${tag2}/" > ${outfile%.*}.tmp4
    mv ${outfile%.*}.tmp4 $outfile
fi

if [ -f ${trait1locifile%.*}.tmp ]; then
    rm ${trait1locifile%.*}.tmp
fi
if [ -f ${trait2locifile%.*}.tmp ]; then
    rm ${trait2locifile%.*}.tmp
fi
rm -f $fdr_clump_loci_file.sorted $fdr_fuma_snp_table.sorted ${outfile%.*}.tmp* ${outfile%.*}_${tag1}_novel_loci.tmp* ${outfile%.*}_${tag2}_novel_loci.tmp*
rm -f ${outfile%.*}_${tag1}_novel_loci.txt ${outfile%.*}_${tag1}_overlap_loci.txt ${outfile%.*}_${tag2}_novel_loci.txt ${outfile%.*}_${tag2}_overlap_loci.txt
