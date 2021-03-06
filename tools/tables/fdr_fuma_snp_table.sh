#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to generate data for supplementary snp table (FUMA+SUMSTAT)
# based on selected snp subset (fdr,r2) from pleioFDR/FUMA analysis WITHOUT
# WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 9 ]; then
  echo "Usage:     sh fdr_fuma_snp_table.sh fdr_clump_snp_file fdr R2 fuma_snp_file sumstat1 sumstat2 tag1 tag2 outfolder"
  echo "Arguments: fdr_clump_snp_file - fdr clumping snp file"
  echo "           fdr - FDR filter for selecting snps"
  echo "           R2 - pleioFDR R2 filter for selecting snps"
  echo "           fuma_snp_file - fuma annotated snp file"
  echo "           sumstat1 - primary summary statistic file"
  echo "           sumstat2 - secondary summary statistic file"
  echo "           tag1 - tag of trait1"
  echo "           tag2 - tag of trait2"
  echo "           outfolder - output folder"
  echo "Example:   sh fdr_fuma_snp_table.sh conj.result.clump.snps.csv 0.1 0.6 snps.txt sumstat/std/trait1.sumstats.gz sumstat/std/trait2.sumstats.gz tag1 tag2 outfolder"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_clump_snp_file=$1
fdr=$2
r2=$3
fuma_snp_file=$4
sumstat1=$5
sumstat2=$6
tag1=$7
tag2=$8
outfolder=$9

mkdir -p $outfolder
outfile=$outfolder/${tag1}_vs_${tag2}_snps.txt

echo "locusnum	CHR	rsID	BP	INDEP_SNP	INDEP_BP	LEAD_SNP	FDR	R2	A1	A2	MAF	gwasP	r2	IndSigSNP	GenomicLocus	nearestGene	dist	func	CADD	RDB	minChrState	commonChrState	posMapFilt	eqtlMapFilt	ciMapFilt	${tag1}_PVAL	${tag1}_Z	${tag1}_BETA	${tag1}_SE	${tag2}_PVAL	${tag2}_Z	${tag2}_BETA	${tag2}_SE	ConcordEffect" > $outfile

cat $fdr_clump_snp_file | awk -v fdr=$fdr -v r2=$r2 'NF==11 && $11<fdr && $7>=r2 {print $6,$1,$2,$5,$4,$3,$8,$11,$7}' | sort -s -k1,1 | uniq > $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt
cat $fuma_snp_file | cut -f2- | sort -s -k1,1 > $outfolder/fuma_snps_${tag1}_${tag2}.txt

sm1=$(basename $sumstat1)
sm2=$(basename $sumstat2)
for sumstat in $sumstat1 $sumstat2; do
    sm=$(basename $sumstat)
    n_z=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n Z | cut -d: -f1`
    if [ `zcat $sumstat | head -n1 | grep OR | wc -l` -gt 0 ]; then
        n_or=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n OR | cut -d: -f1`
        if [ `zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep SE | grep -v NCASE | wc -l` -gt 0 ]; then
            n_se=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n SE | grep -v NCASE | cut -d: -f1`
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_or=$n_or -v n_se=$n_se '{print $1,$5,$6,$4,$n_z,log($n_or),$n_se/$n_or}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        else
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_or=$n_or '{if($n_z==0) print $1,$5,$6,$4,$n_z,log($n_or),0; else print $1,$5,$6,$4,$n_z,log($n_or),log($n_or)/$n_z}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        fi
    elif [ `zcat $sumstat | head -n1 | grep BETA | wc -l` -gt 0 ]; then
        n_beta=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n BETA | cut -d: -f1`
        if [ `zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep SE | grep -v NCASE | wc -l` -gt 0 ]; then
            n_se=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n SE | grep -v NCASE | cut -d: -f1`
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_beta=$n_beta -v n_se=$n_se '{print $1,$5,$6,$4,$n_z,$n_beta,$n_se}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        else
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_beta=$n_beta '{if($n_z==0) print $1,$5,$6,$4,$n_z,$n_beta,0; else print $1,$5,$6,$4,$n_z,$n_beta,$n_beta/n_z}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        fi
    elif [ `zcat $sumstat | head -n1 | grep LOGODDS | wc -l` -gt 0 ]; then
        n_beta=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n LOGODDS | cut -d: -f1`
        if [ `zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep SE | grep -v NCASE | wc -l` -gt 0 ]; then
            n_se=`zcat $sumstat | head -n1 | sed 's/	/\n/g' | grep -n SE | grep -v NCASE | cut -d: -f1`
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_beta=$n_beta -v n_se=$n_se '{print $1,$5,$6,$4,$n_z,$n_beta,$n_se}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        else
            zcat $sumstat | tail -n +2 | awk -v n_z=$n_z -v n_beta=$n_beta '{if($n_z==0) print $1,$5,$6,$4,$n_z,$n_beta,0; else print $1,$5,$6,$4,$n_z,$n_beta,$n_beta/n_z}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
        fi
    else
        zcat $sumstat | tail -n +2 | awk -v n_z=$n_z '{print $1,$5,$6,$4,$n_z,"NA NA"}' | sort -s -k1,1 > $outfolder/${sm%%.*}.txt
    fi
done

#join -1 1 -2 1 -a 1 $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt $outfolder/fuma_snps_${tag1}_${tag2}.txt | awk '{if(NF==9) print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" - - - - - - - - - - - - - - - - - - -"; else print $0}' | cut -d' ' -f1-9,12- > $outfolder/fuma_snps_${tag1}_${tag2}.tmp
join -1 1 -2 1 $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt $outfolder/fuma_snps_${tag1}_${tag2}.txt | cut -d' ' -f1-9,12- > $outfolder/fuma_snps_${tag1}_${tag2}.tmp

join -1 1 -2 1 -a 1 $outfolder/fuma_snps_${tag1}_${tag2}.tmp $outfolder/${sm1%%.*}.txt > $outfolder/fuma_snps_${tag1}_${tag2}.tmp2
join -1 1 -2 1 -a 1 $outfolder/fuma_snps_${tag1}_${tag2}.tmp2 $outfolder/${sm2%%.*}.txt | sort -k1,1 -u | awk '{if($10=="-" || $11=="-") print $2,$3,$1,$4,$5,$6,$7,$8,$9,$27,$28,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38; else print $2,$3,$1,$4,$5,$6,$7,$8,$9,$11,$10,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38}'| awk '{if($10==$28 && $11==$27&& $31!="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,-$30,-$31,$32,$33,$34,$35,$36,$37,$38; else if($10==$28 && $11==$27 && $31=="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,-$30,"NA",$32,$33,$34,$35,$36,$37,$38; else print $0}' | awk '{if($10==$34 && $11==$33 && $37!="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,-$36,-$37,$38; else if($10==$34 && $11==$33 && $37=="NA") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,-$36,"NA",$38; else print $0}' | awk '{if($30<0 && $36<0 || $30>0 && $36>0) print $0,"TRUE"; else print $0,"FALSE"}' | cut -d' ' -f1-26,29-32,35- | sed 's/ /	/g' | sort -n -k2,2 -k4,4 >> $outfile

cut -f1-12,14- $outfile > $outfile.tmp
mv $outfile.tmp $outfile

n_loci=`tail -n +2 $outfile | cut -f1 | sort | uniq | wc -l`
n_chr=`tail -n +2 $outfile | cut -f2 | sort | uniq | wc -l`
n_rsid=`tail -n +2 $outfile | cut -f3 | wc -l`
n_indep=`tail -n +2 $outfile | cut -f5 | sort -s | uniq | wc -l`
n_near=`tail -n +2 $outfile | cut -f16 | sort -s | uniq | wc -l`
n_concord=`cut -f34 $outfile | grep TRUE | wc -l`
echo $n_loci $n_chr $n_rsid $n_indep $n_near $n_concord | awk '{print $1,$2,$3,$4,$5,$6/$3}' OFMT="%.4f" | awk '{print $1,$2,$3,$4,$5,$6*100"%"}' | awk '{print $1"\t"$2"\t"$3"\t\t"$4"\t\t\t\t\t\t\t\t\t\t\t"$5"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"$6}' >> $outfile

rm -f $outfolder/fdr_clump_snps_${tag1}_${tag2}.txt $outfolder/fuma_snps_${tag1}_${tag2}.txt $outfolder/${sm1%%.*}.txt $outfolder/${sm2%%.*}.txt $outfolder/fuma_snps_${tag1}_${tag2}.tmp $outfolder/fuma_snps_${tag1}_${tag2}.tmp2
