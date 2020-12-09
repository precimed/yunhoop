#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary snp table (FUMA+SUMSTAT)
# based on selected snp subset (fdr,r2) from cFDR/FUMA analysis, but WITHOUT
# ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 7 ]; then
  echo "Usage:     sh fdr_fuma_snp_table.sh cfdr_clump_snp_file fdr r2 fuma_snp_file sumstat1 sumstat2 outfile"
  echo "Arguments: cfdr_clump_snp_file - file that contains cfdr clumping snp info"
  echo "           fdr - FDR filter for selecting snps"
  echo "           r2 - r2 filter for selecting snps"
  echo "           fuma_snp_file - FUMA annotated snp file"
  echo "           sumstat1 - primary summary statistic file"
  echo "           sumstat2 - secondary summary statistic file"
  echo "           outfile - output file"
  echo "Example:   sh fdr_fuma_snp_table.sh PGC_BIP_2016_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.snps.csv 0.1 0.6 snp2gene/BIP_vs_MOOD_conj_005/snps.txt sumstat/std/PGC_BIP_2016.sumstats.gz sumstat/std/UKB_MOOD_2019.sumstats.gz snp2gene/BIP_vs_MOOD_snps.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

cfdr_clump_snp_file=$1
fdr=$2
r2=$3
fuma_snp_file=$4
sumstat1=$5
sumstat2=$6
outfile=$7
outfolder=`dirname $outfile`

tag1=`echo $sumstat1 | cut -d'_' -f2`
tag2=`echo $sumstat2 | cut -d'_' -f2`

echo "rsID	chr	pos	LEAD_SNP	FDR	non_effect_allele	effect_allele	r2	IndSigSNP	GenomicLocus	nearestGene	dist	func	CADD	RDB	minChrState	commonChrState	posMapFilt	eqtlMapFilt	ciMapFilt	${tag1}_PVAL	${tag1}_Z	${tag1}_BETA	${tag2}_PVAL	${tag2}_Z	${tag2}_BETA	Concordant_Effect" > $outfile

cat $cfdr_clump_snp_file | awk -v fdr=$fdr -v r2=$r2 'NF==11 && $11<fdr && $7>=r2 {print $6,$8,$11}' | sort -s -k1,1 | uniq > $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt

cat $fuma_snp_file | cut -f2-6,9- | sort -s -k1,1 > $outfolder/fuma_snps_${tag1}_${tag2}.txt

sm1=$(basename $sumstat1)
sm2=$(basename $sumstat2)
for sumstat in $sumstat1 $sumstat2; do
    sm=$(basename $sumstat)
    n_z=`zcat $sumstat | head -n1 | sed 's/\t/\n/g' | grep -n Z | cut -d: -f1`
    if [ `zcat $sumstat | head -n1 | grep OR | wc -l` -gt 0 ]; then
        n_or=`zcat $sumstat | head -n1 | sed 's/\t/\n/g' | grep -n OR | cut -d: -f1`
        zcat $sumstat | awk -v n_z=$n_z -v n_or=$n_or '{print $1,$5,$6,$4,$n_z,log($n_or)}' | sort -s -k1,1 > ${sm%%.*}.txt
    elif [ `zcat $sumstat | head -n1 | grep BETA | wc -l` -gt 0 ]; then
        n_beta=`zcat $sumstat | head -n1 | sed 's/\t/\n/g' | grep -n BETA | cut -d: -f1`
        zcat $sumstat | awk -v n_z=$n_z -v n_beta=$n_beta '{print $1,$5,$6,$4,$n_z,$n_beta}' | sort -s -k1,1 > ${sm%%.*}.txt
    else
        zcat $sumstat | awk -v n_z=$n_z '{print $1,$5,$6,$4,$n_z,"NA"}' | sort -s -k1,1 > ${sm%%.*}.txt
    fi
done

join -1 1 -2 1 $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt $outfolder/fuma_snps_${tag1}_${tag2}.txt > $outfolder/fuma_snps_${tag1}_${tag2}.tmp
join -1 1 -2 1 $outfolder/fuma_snps_${tag1}_${tag2}.tmp ${sm1%%.*}.txt > $outfolder/fuma_snps_${tag1}_${tag2}.tmp2
join -1 1 -2 1 $outfolder/fuma_snps_${tag1}_${tag2}.tmp2 ${sm2%%.*}.txt | awk '{if ($6==$21 && $7==$22 && $25!="NA") print $1,$4,$5,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,-$24,-$25,$26,$27,$28,$29,$30; else if ($6==$21 && $7==$22 && $25=="NA") print $1,$4,$5,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,-$24,"NA",$26,$27,$28,$29,$30; else print $1,$4,$5,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' | awk '{if ($6==$26 && $7==$27 && $30!="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,-$29,-$30; else if ($6==$26 && $7==$27 && $30=="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,-$29,"NA"; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30}' | awk '{if ($25<0 && $30<0 || $25>0 && $30>0) print $0,"TRUE"; else print $0,"FALSE"}' | cut -d' ' -f1-20,23-25,28- | sed 's/ / /g' | sort -n -k2,2 -k3,3 >> $outfile

rm $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt $outfolder/fuma_snps_${tag1}_${tag2}.txt ${sm1%%.*}.txt ${sm2%%.*}.txt $outfolder/fuma_snps_${tag1}_${tag2}.tmp $outfolder/fuma_snps_${tag1}_${tag2}.tmp2
