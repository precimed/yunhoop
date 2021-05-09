#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to ensemble NHGRI-EBI GWAS Catalog to make a clean file
# for novelty checking of loci WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 3 ]; then
  echo "Usage: sh compile_gwasc.sh gwasc dbsnp_map_folder reffile outfile"
  echo "Arguments: gwasc - file of gwascatalog"
  echo "           dbsnp_map_folder - top folder to hold built dbsnp maps within json and vcf subfolders"
  echo "           outfile - output file"
  echo "Example: sh compile_gwasc.sh gwas_catalog_v1.0-associations_e100_r2021-04-12.tsv $lf/ncbi/b153 $data/ref/9545380.ref gwas_catalog_v1.0-associations_e100_r2021-04-12.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
dbsnp_map_folder=$2
outfile=$3
chainfile_hg1538=$yc/software/liftover/ucsc/hg15ToHg38.over.chain.gz
chainfile_hg1619=$yc/software/liftover/ucsc/hg16ToHg19.over.chain.gz
chainfile_hg1638=$yc/software/liftover/ucsc/hg16ToHg38.over.chain.gz
chainfile_hg1719=$yc/software/liftover/ucsc/hg17ToHg19.over.chain.gz
chainfile_hg1738=$yc/software/liftover/ucsc/hg17ToHg38.over.chain.gz
chainfile_hg1819=$yc/software/liftover/ucsc/hg18ToHg19.over.chain.gz
chainfile_hg1838=$yc/software/liftover/ucsc/hg18ToHg38.over.chain.gz
chainfile_hg1938=$yc/software/liftover/ucsc/hg19ToHg38.over.chain.gz
chainfile_hg3819=$yc/software/liftover/ucsc/hg38ToHg19.over.chain.gz
liftover=$lf/liftOver
outfolder=$(dirname $outfile)
reffile=$data/ref/9545380.ref

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
cut -f7-8,12-13,21-24 $gwasc | awk -F '\t' '{print $3,$4,$6,$5,$7,$8,$1,$2}' OFS='\t' | tr -dc '\0-\177' | sed 's/&beta;/Beta /g' > $gwasc.txt

#filter snps with irregular ids
echo "filter irregular snps with empty chr"
cut -f1-3 $gwasc.txt | awk -F '\t' '$1==""' | grep -v ":" | grep -i -v rs | grep -i -v chr | cut -f3 > $gwasc.tmp
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $gwasc.tmp $gwasc.txt > $gwasc.tmp2
mv $gwasc.tmp2 $gwasc.txt

if [ `cut -f3 $gwasc.txt | grep rs | grep : | wc -l` -gt 0 ]; then
   echo "NOTE: IDs combined with rs and without rs"
   cut -f3 $gwasc.txt | grep rs | grep :
   exit 1
fi

#clean snp ids
echo "clean snps"
tail -n +2 $gwasc.txt | cut -f1-5 | sed 's/\./:/g' | sed 's/\s+//g' | sed 's/ ; /;/g' | sed 's/; /;/g' | sed 's/ ;/;/g' | sed 's/ , /;/g' | sed 's/, /;/g' | sed 's/ ,/;/g' | sed 's/,/;/g' | sed 's/ \/ /,/g' | sed 's/ \//,/g' | sed 's/\/ /,/g' | sed 's/\//,/g' | sed 's/ x /x/gi' > $gwasc.tmp1
tail -n +2 $gwasc.txt | cut -f6 | cut -d'_' -f1 | cut -d'-' -f1 | sed 's/[a-zA-Z]$//' | sed 's/^/rs/' | sed 's/^rs$/-/' > $gwasc.tmp2
tail -n +2 $gwasc.txt | cut -f7-8 > $gwasc.tmp3
for((i=0; i<=22; i++)); do
    sed "s/chr${i}: /chr${i}:/gi" $gwasc.tmp1 > $gwasc.tmp10
    sed "s/chr${i}_/chr${i}:/gi" $gwasc.tmp10 > $gwasc.tmp1
done
cut -f1-3 $gwasc.tmp1 | sed 's/chr://gi' | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/		.*psy_/		/gi' | sed 's/hg18_/hg18/gi' > $gwasc.tmp10
cut -f4-5 $gwasc.tmp1 | sed 's/chr://gi' | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/^.*psy_/psy_/gi' > $gwasc.tmp11
paste -d'	' $gwasc.tmp10 $gwasc.tmp11 $gwasc.tmp2 $gwasc.tmp3 > $gwasc.tmp
rm -f $gwasc.tmp1* $gwasc.tmp2 $gwasc.tmp3

if [ `awk -F '\t' '($3~/x/ || $3~/;/)&&$6!="-" {print $3,$6}' OFS='\t' $gwasc.tmp | wc -l` -gt 0 ]; then
    echo "multi ids with single current id"
    awk -F '\t' '($3~/x/ || $3~/;/)&&$6!="-" {print $3,$6}' OFS='\t' $gwasc.tmp
    exit 1
fi

#process records with only single or multiple ids (chr:pos or rs)
echo "process records with only single or multiple ids"
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/' > $gwasc.tmp1
rm -f $gwasc.tmp12
cut -f3- $gwasc.tmp1 | while read line; do
    rm -f $gwasc.tmp10
    rest=`echo "$line" | cut -f3-`
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | while read id; do
        snp=`echo $id | sed 's/del-//gi' | cut -d'_' -f1 | cut -d'-' -f1 | cut -d':' -f1-2 | sed 's/[a-zA-Z]$//'`
        if [ `echo $snp | grep ':' | wc -l` -gt 0 ]; then
            chr=`echo $snp | cut -d':' -f1`
            pos=`echo $snp | cut -d':' -f2`
            echo $chr'	'$pos'	'$snp'	'$id'	'"$rest" >> $gwasc.tmp10
        else
            echo '-	-	'$snp'	'$id'	'"$rest" >> $gwasc.tmp10
        fi
    done 
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp11
    n0=`cat $gwasc.tmp10 | wc -l`
    n1=`cat $gwasc.tmp11 | wc -l`
    if [ $n0 -ne $n1 ]; then
        echo "NOTE: snp number ($n0) does not match risk-allele id number ($n1)"
        echo $line
        exit 1
    fi
    paste -d'	' $gwasc.tmp10 $gwasc.tmp11 | awk -F '\t' '{print $1,$2,$3,$4,$9,$5,$6,$7,$8}' OFS='\t' >> $gwasc.tmp12
    rm -f $gwasc.tmp10 $gwasc.tmp11
done
mv $gwasc.tmp12 $gwasc.tmp1

#process records as single snps with chr:pos and id
echo "process records as single snp with chr:pos and id (chr:pos or rs)"
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/ && $3!~/x/' $gwasc.tmp > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/del-//gi'| sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f3- $gwasc.tmp2 > $gwasc.tmp23
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 > $gwasc.tmp2
rm -f $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#process records associated with multiple snps with chr:pos and ids
echo "process records with multiple snps with chr:pos and ids"
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/ || $3~/x/)' $gwasc.tmp > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp31
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp32
    echo "$line" | cut -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
    echo "$line" | cut -f4 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | awk '$1~/:/ || $1~/rs/ || ($1!~/rs/ && $1!~/[a-zA-Z]/)' > $gwasc.tmp34
    rest=`echo "$line" | cut -f5-`
    n1=`cat $gwasc.tmp31 | wc -l`
    n2=`cat $gwasc.tmp32 | wc -l`
    n3=`cat $gwasc.tmp33 | wc -l`
    n4=`cat $gwasc.tmp34 | wc -l`
    if [ $n1 -ne $n2 ]; then
        echo $line
        echo "NOTE: chr number ($n1) does not match pos number ($n2)"
        exit 1
    fi
    if [ $n3 -ne $n4 ]; then
        echo $line
        echo "NOTE: snp number ($n3) does not match risk-allele id number ($n4)"
        exit 1
    fi
    if [ $n1 -gt $n3 ]; then
        echo $line
        echo "NOTE: chr number ($n1) is greater than snp number ($n3)"
        exit 1
    fi
    # a) number of chr:pos matches number of ids
    if [ $n3 -eq $n1 ] && [ `echo "$line" | cut -f3 | grep 'x' | wc -l` -gt 0 ]; then
        rm -f $gwasc.tmp35
        for ((i=1; i<=$n3; i++)); do
            echo "$rest" >> $gwasc.tmp35
        done
    # b) single chr:pos with multiple merging rs
    elif [ $n1 -eq 1 ] && [ $n3 -gt $n1 ] && [ `echo "$line" | cut -f3 | grep ',' | wc -l` -gt 0 ]; then
        for((i=2; i<=$n3; i++)); do
            head -n1 $gwasc.tmp31 >> $gwasc.tmp31
            head -n1 $gwasc.tmp32 >> $gwasc.tmp32
        done
        rm -f $gwasc.tmp35
        for ((i=1; i<=$n3; i++)); do
            echo "$rest" >> $gwasc.tmp35
        done
    else
    #    paste -d':' $gwasc.tmp31 $gwasc.tmp32 > $gwasc.tmp34
    #    for((i=1; i<=$n3; i++)); do
    #        echo '-' >> $gwasc.tmp31
    #        echo '-' >> $gwasc.tmp32
    #    done
    #    cat $gwasc.tmp33 >> $gwasc.tmp34
    #    mv $gwasc.tmp34 $gwasc.tmp33
    #    rm -f $gwasc.tmp34
    #    for ((i=1; i<=$((n1+n3)); i++)); do
    #        echo "$study" >> $gwasc.tmp34
    #    done

    # c) number of chr:pos doesn't match number of ids (take as ids with empty chr:pos)
        rm -f $gwasc.tmp31
        rm -f $gwasc.tmp32
        rm -f $gwasc.tmp35
        for((i=1; i<=$n3; i++)); do
            echo '-' >> $gwasc.tmp31
            echo '-' >> $gwasc.tmp32
            echo "$rest" >> $gwasc.tmp35
        done
    fi
    sed 's/del-//gi' $gwasc.tmp33 | sed 's/[a-zA-Z]$//' > $gwasc.tmp36
    if [ `grep : $gwasc.tmp36 | wc -l` -gt 0 ]; then
        echo "NOTE: snp without rsid is found"
        echo $line
        exit 1
    fi
    paste -d'	' $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp36 $gwasc.tmp33 $gwasc.tmp34 $gwasc.tmp35 >> $gwasc.tmp30
    rm -f $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp36 $gwasc.tmp33 $gwasc.tmp34 $gwasc.tmp35
done
mv $gwasc.tmp30 $gwasc.tmp3

cat $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3 | sed 's/^23	/X	/' | sed 's/	23:/	X:/' | sed 's/^24	/Y	/' | sed 's/	24:/	Y:/' | sort | uniq > $outfile
rm -f $gwasc.tmp $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3

#update chr:pos of records with only rs in light of available snps
echo "update chr:pos of records with only rs in light of available snps"
grep -v ^- $outfile > $outfile.tmp
cut -f1-3 $outfile.tmp | sort | uniq | sort -s -k3,3 > $outfile.tmp1
grep ^- $outfile | sort | uniq | sort -s -k3,3 > $outfile.tmp2
join -1 3 -2 3 -a 1 -t '	' $outfile.tmp2 $outfile.tmp1 | awk -F '\t' '{if (NF>9) print $10,$11,$1,$4,$5,$6,$7,$8,$9; else print $2,$3,$1,$4,$5,$6,$7,$8,$9}' OFS='\t' >> $outfile.tmp

#create initial set (OUTPUT: CHR POS SNP ID RISK_ALLELE MERGED CURRENT_ID STUDY TRAIT FLAG)
echo "compile initial set"
grep -v ^- $outfile.tmp | awk -F '\t' '$1~/^[0-9]+$/' | sort -n -k1,1 -k2,2 > $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1!~/^[0-9]+$/' | sort -k1,1 >> $outfile.tmp2
grep ^- $outfile.tmp >> $outfile.tmp2
echo 'CHR	POS	SNP	ID	RISK_ALLELE	MERGED	CURRENT_ID	STUDY	TRAIT	FLAG' > $outfile
#flag: 0-snp with only rs; 1-snp with chr:pos; 2-snp with both rs and chr:pos
awk -F '\t' '{if($1=="-") print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0"; else if($1!="-"&&$3~/:/) print $1,$2,$3,$4,$5,$6,$7,$8,$9,"1"; else if($1!="-"&&$3~/rs/) print $1,$2,$3,$4,$5,$6,$7,$8,$9,"2"; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,"x"}' OFS='\t' $outfile.tmp2 >> $outfile

if [ `cut -f1-3 $outfile | grep : | sed 's/hg18//g' | awk '$1":"$2!=$3' | wc -l` -gt 0 ]; then
    echo "NOTE: SNPs with inconsistent chr:pos" 
    cut -f1-3 $outfile | grep : | sed 's/hg18//g' | awk '$1":"$2!=$3'
    exit 1
fi
if [ `cut -f1-3 $outfile | grep '-' | grep ':' | wc -l` -gt 0 ]; then
    echo "NOTE: SNPs with chr:pos to split" 
    cut -f1-3 $outfile | grep '-' | grep ':'
    exit 1
fi

#merge SNP#3 to CURRENT_ID#7
head -n1 $outfile > ${outfile%.*}_cur.csv
tail -n +2 $outfile | awk -F '\t' '{if($7=="-") print $1,$2,$3,$4,$5,$6,$3,$8,$9,$10; else print $0}' OFS='\t' | sort | uniq >> ${outfile%.*}_cur.csv
tail -n +2 ${outfile%.*}_cur.csv | cut -f1,2,3,7 | sed 's/[0-9]//g' | sed 's/rs//g' | sed 's/hg//g' | sed 's/[X-Y]//g' | sed 's/://g' | sed 's/-//g' | sort | uniq
cut -f7 ${outfile%.*}_cur.csv | grep rs | sort | uniq > ${outfile%.*}_rs.txt
rm -f $outfile.tmp*
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build dbSNP misc maps
if [ ! -f $dbsnp_map_folder/json/dbsnp_merge_map.txt ]; then
    sh $(dirname $0)/../dbsnp/build_dbsnp_maps.sh $dbsnp_map_folder N N Y
    cat $dbsnp_map_folder/json/dbsnp_unsupported.txt $dbsnp_map_folder/json/dbsnp_withdrawn.txt $dbsnp_map_folder/json/dbsnp_nosnppos.txt > $dbsnp_map_folder/json/dbsnp_unavailable.txt
fi

#build dbSNP chr maps
if [ ! -f $dbsnp_map_folder/vcf/dbsnp_b38_chr1.txt ]; then
    sh $(dirname $0)/../dbsnp/build_dbsnp_maps.sh $dbsnp_map_folder N Y N
fi
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build ref snp merge map (based on ref rs list [9545380.ref] and dbsnp_merge_map.txt)
#if [ ! -f ${reffile%.*}_snp_merge_map.txt ]; then
#    tail -n +2 $reffile | cut -f2 > ${reffile%.*}_rs.txt
#    sh $(dirname $0)/build_snp_merge_map.sh $dbsnp_map_folder/json/dbsnp_merge_map.txt ${reffile%.*}_rs.txt ${reffile%.*}_snp_merge_map.txt
#    sort -s -k2,2 $reffile > ${reffile%.*}.sorted
#    sort -s -k1,1 ${reffile%.*}_snp_merge_map.txt > ${reffile%.*}_snp_merge_map.sorted
#    echo 'CHR	CURRENT_ID	SNP	BP	A1	A2	MERGED' > ${reffile%.*}.txt
#    join -1 2 -2 1 -a 1 ${reffile%.*}.sorted ${reffile%.*}_snp_merge_map.sorted -t '	' | grep -v SNP | awk -F '\t' '{if(NF>6) print $2,$7,$1,$4,$5,$6,$8; else print $2,$1,$1,$4,$5,$6,"0"}' OFS='\t' | sort -n -k1,1 -k4,4 >> ${reffile%.*}.txt
#    rm -f ${reffile%.*}.sorted ${reffile%.*}_snp_merge_map.sorted
#fi

#get updated dbsnp b37/b38 chr:pos after merging
#echo 'CURRENT_RS	CHR	SNP	BP	A1	A2	MERGED	POS	REF	ALT	ID	CURRENT_ID' > ${reffile%.*}_b37.txt
#echo 'CURRENT_RS	CHR	SNP	BP	A1	A2	MERGED	POS	REF	ALT	ID	CURRENT_ID' > ${reffile%.*}_b38.txt
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
#    echo $i
#    grep "^$i	" ${reffile%.*}.txt | sort -s -k2,2 > ${reffile%.*}_chr$i.tmp
#    awk -F '\t' 'NR==FNR {D[$2]++; next} ($1 in D)' ${reffile%.*}_chr$i.tmp $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt | sort -s -k1,1 > $(dirname $reffile)/dbsnp_b37_chr$i.tmp
#    awk -F '\t' 'NR==FNR {D[$2]++; next} ($1 in D)' ${reffile%.*}_chr$i.tmp $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt | sort -s -k1,1 > $(dirname $reffile)/dbsnp_b38_chr$i.tmp
#    join -1 2 -2 1 ${reffile%.*}_chr$i.tmp $(dirname $reffile)/dbsnp_b37_chr$i.tmp -t '	' | cut -f1-7,10-12 | awk -F '\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2":"$4,$2":"$8}' OFS='\t' | sort -n -k2,2 -k4,4 >> ${reffile%.*}_b37.txt
#    join -1 2 -2 1 ${reffile%.*}_chr$i.tmp $(dirname $reffile)/dbsnp_b38_chr$i.tmp -t '	' | cut -f1-7,10-12 | awk -F '\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2":"$4,$2":"$8}' OFS='\t' | sort -n -k2,2 -k4,4 >> ${reffile%.*}_b38.txt
#done
#rm -f ${reffile%.*}_chr*.tmp $(dirname $reffile)/dbsnp_b*_chr*.tmp

#build gwasc snp merge map (based on gwasc rs list and dbsnp_merge_map.txt)
sh $(dirname $0)/build_snp_merge_map.sh $dbsnp_map_folder/json/dbsnp_merge_map.txt ${outfile%.*}_rs.txt $outfolder/gwasc_snp_merge_map.txt

#update gwasc rs list with merged-into snps (based on gwasc_snp_merge_map.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_rs.txt > ${outfile%.*}_rs_merged.txt
cut -f2 $outfolder/gwasc_snp_merge_map.txt | sort | uniq >> ${outfile%.*}_rs_merged.txt

#update gwasc rs list with available snps (based on dbsnp_unavailable.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_rs_merged.txt | sort | uniq > ${outfile%.*}_rs_avail.txt

#update in light of gwasc snp merge map
awk -F '\t' 'NR==FNR {D[$1]++; next} !($7 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_cur.csv | awk -F '\t' '{print $1,$2,$7,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' > ${outfile%.*}_merged.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($7 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_cur.csv | sort -s -k7,7 > ${outfile%.*}.tmp
sort -s -k1,1 $outfolder/gwasc_snp_merge_map.txt > $outfolder/gwasc_snp_merge_map.tmp
join -1 7 -2 1 -a 1 -t '	' ${outfile%.*}.tmp $outfolder/gwasc_snp_merge_map.tmp | awk -F '\t' '{if(NF>10) print $2,$3,$11,$4,$5,$6,$7,$1,$8,$9,$10; else print $2,$3,$1,$4,$5,$6,$7,$1,$8,$9,$10}' OFS='\t' >> ${outfile%.*}_merged.txt
rm -f ${outfile%.*}.tmp $outfolder/gwasc_snp_merge_map.tmp

#update in light of dbsnp_unavailable.txt
head -n1 ${outfile%.*}_merged.txt > ${outfile%.*}_avail.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_merged.txt | tail -n +2 | sort | uniq >> ${outfile%.*}_avail.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($3 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_merged.txt | cut -f3 | sort | uniq > ${outfile%.*}_unavail.txt

#snps with only chr:pos
cut -f3 ${outfile%.*}_avail.txt | grep : | grep -v hg |  sort -s | uniq > ${outfile%.*}_bp.txt

if [ `cut -f3,8 ${outfile%.*}_avail.txt | sort | uniq | cut -f2 | sort | uniq -c | awk '$1>1' | wc -l` -gt 0 ]; then
    echo "snps merged more than once:"
    cut -f3,8 ${outfile%.*}_avail.txt | sort | uniq | cut -f2 | sort | uniq -c | awk '$1>1'
fi
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build gwasc snp chr map (based on available rs list and dbsnp chr map)
#sh $(dirname $0)/build_snp_map.sh $dbsnp_map_folder ${outfile%.*}_rs_avail.txt $outfolder

#build b38 and b37 map for snps with rs
#rm -f $outfolder/gwasc_A_b38_dbsnp_rs_map.txt
#rm -f $outfolder/gwasc_A_b37_dbsnp_rs_map.txt
rm -f $outfolder/gwasc_B_b38_dbsnp_rs_map.txt
rm -f $outfolder/gwasc_B_b37_dbsnp_rs_map.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    echo $i
    #awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b38.txt >> $outfolder/gwasc_A_b38_dbsnp_rs_map.txt
    #awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b37.txt >> $outfolder/gwasc_A_b37_dbsnp_rs_map.txt
    #input: $2:0-based coordinate $9:1-based coordinate
    #output: $2:1-based coordinate $9:0-based coordinate
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_B_b38_dbsnp_rs_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_B_b37_dbsnp_rs_map.txt
done

awk -F '\t' '{print $2"_"$1,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' $outfolder/gwasc_B_b38_dbsnp_rs_map.txt > $outfolder/gwasc_B_b38_dbsnp_rs_map.tmp
awk -F '\t' '{print $3"_"$1,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' $outfolder/gwasc_B_b37_dbsnp_rs_map.txt > $outfolder/gwasc_B_b37_dbsnp_rs_map.tmp
cut -f1-3 ${outfile%.*}_avail.txt | grep rs | grep '-' | cut -f3 | sort | uniq > ${outfile%.*}_snp_avail_0.txt
cut -f1-3 ${outfile%.*}_avail.txt | grep rs | grep -v '-' | awk -F '\t' '{print $1":"$2"_"$3,$3}' OFS='\t' | sort | uniq > ${outfile%.*}_snp_avail_2.txt
cut -f1-3 ${outfile%.*}_avail.txt | grep rs | grep -v '-' | awk -F '\t' '{print $1"_"$3,$3}' OFS='\t' | sort | uniq > ${outfile%.*}_snp_avail_3.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($1 in D)' ${outfile%.*}_snp_avail_3.txt ${outfile%.*}_snp_avail_0.txt > ${outfile%.*}_snp_avail_1.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_snp_avail_1.txt $outfolder/gwasc_B_b38_dbsnp_rs_map.tmp | cut -f2- > $outfolder/gwasc_B_b38_dbsnp_rs_map2.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_snp_avail_2.txt $outfolder/gwasc_B_b38_dbsnp_rs_map.tmp | cut -f2- >> $outfolder/gwasc_B_b38_dbsnp_rs_map2.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_snp_avail_1.txt $outfolder/gwasc_B_b37_dbsnp_rs_map.tmp | cut -f2- > $outfolder/gwasc_B_b37_dbsnp_rs_map2.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_snp_avail_3.txt $outfolder/gwasc_B_b37_dbsnp_rs_map.tmp | cut -f2- >> $outfolder/gwasc_B_b37_dbsnp_rs_map2.txt
rm -f $outfolder/gwasc_*_b3*_dbsnp_rs_map.tmp
rm -f ${outfile%.*}_snp_avail_*.txt
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build b38 and b37 map for snps with only chr:pos
rm -f $outfolder/gwasc_B_b38_dbsnp_chr_map.txt
rm -f $outfolder/gwasc_B_b37_dbsnp_chr_map.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    echo $i
    #0-based coordinate
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr${i}.txt >> $outfolder/gwasc_B_b38_dbsnp_chr_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr${i}.txt >> $outfolder/gwasc_B_b37_dbsnp_chr_map.txt
    #awk -F '\t' 'NR==FNR {D[$1]++; next} ($9 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_B_b38_dbsnp_chr_map.txt
    #awk -F '\t' 'NR==FNR {D[$1]++; next} ($9 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_B_b37_dbsnp_chr_map.txt
done
fi

tag='B'
run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#chr:pos shared between b38 and b37 are taken as b38 (since it currently mapped to genome assembly GRCh38)
awk -F '\t' 'NR==FNR {D[$2]++; next} ($2 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map.txt $outfolder/gwasc_${tag}_b38_dbsnp_chr_map.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map1.txt
cut -f2 $outfolder/gwasc_${tag}_b38_dbsnp_chr_map1.txt | sort | uniq > $outfolder/gwasc_${tag}_b3738_dbsnp_comm.txt

#chr:pos only seen with b38 taken as b38
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map.txt $outfolder/gwasc_${tag}_b38_dbsnp_chr_map.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map2.txt

#chr:pos only seen with b37 taken as b37
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map.txt $outfolder/gwasc_${tag}_b37_dbsnp_chr_map.txt > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map2.txt

cat $outfolder/gwasc_${tag}_b38_dbsnp_chr_map1.txt $outfolder/gwasc_${tag}_b38_dbsnp_chr_map2.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt
cat $outfolder/gwasc_${tag}_b37_dbsnp_chr_map2.txt > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt
rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.tmp
rm -f $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.tmp
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    echo $i
    #0-based coordinate
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr${i}.txt >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.tmp
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr${i}.txt >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.tmp
done
cat $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.tmp >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt
cat $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.tmp >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt
sort -n -k3,3 -k4,4 $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.tmp
sort -n -k3,3 -k4,4 $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.tmp
mv $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.tmp $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt
mv $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.tmp $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#chr:pos not seen in both (check for hg18)
cat $outfolder/gwasc_${tag}_b37_dbsnp_chr_map.txt $outfolder/gwasc_${tag}_b38_dbsnp_chr_map.txt | cut -f2 | sort | uniq > $outfolder/gwasc_${tag}_b3738_dbsnp_summ.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $outfolder/gwasc_${tag}_b3738_dbsnp_summ.txt ${outfile%.*}_bp.txt | sort | uniq > $outfolder/gwasc_${tag}_b3738_dbsnp_miss.txt
cat $outfolder/gwasc_${tag}_b3738_dbsnp_miss.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg18.txt
cut -f3 ${outfile%.*}_avail.txt | grep hg18 | sed 's/hg18//' | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),"hg18"$1":"$2}' OFS='\t' >> ${outfile%.*}_hg18.txt
if [ -s ${outfile%.*}_hg18.txt ]; then
    #liftover from hg18 to hg19
    $liftover ${outfile%.*}_hg18.txt $chainfile_hg1819 ${outfile%.*}_lifted_hg18_to_hg19.bed ${outfile%.*}_unlifted_hg18_to_hg19.bed
    rm -f $outfolder/gwasc_${tag}_b37_dbsnp_chr_map3.txt
    cut -f1-2 ${outfile%.*}_lifted_hg18_to_hg19.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map3.txt
        #awk -v snp=$snp -F '\t' '$9==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map3.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg18_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg18_hg19.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map3.txt ${outfile%.*}_hg18_hg19.txt | cut -f1 > ${outfile%.*}_hg18_unmapped_b37.txt

    #liftover from hg18 to hg38
    $liftover ${outfile%.*}_hg18.txt $chainfile_hg1838 ${outfile%.*}_lifted_hg18_to_hg38.bed ${outfile%.*}_unlifted_hg18_to_hg38.bed
    rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map3.txt
    cut -f1-2 ${outfile%.*}_lifted_hg18_to_hg38.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map3.txt
        #awk -v snp=$snp -F '\t' '$9==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt | awk -F '\t' '{print $1,$9,$3,$4,$5,$6,$7,$8,$2}' OFS='\t' >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map3.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg18_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg18_hg38.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map3.txt ${outfile%.*}_hg18_hg38.txt | cut -f1 > ${outfile%.*}_hg18_unmapped_b38.txt
fi

#chr:pos neither seen in both (check for hg17)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg18_unmapped_b38.txt ${outfile%.*}_hg18_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg17.txt
if [ -s ${outfile%.*}_hg17.txt ]; then
    #liftover from hg17 to hg19
    $liftover ${outfile%.*}_hg17.txt $chainfile_hg1719 ${outfile%.*}_lifted_hg17_to_hg19.bed ${outfile%.*}_unlifted_hg17_to_hg19.bed
    rm -f $outfolder/gwasc_${tag}_b37_dbsnp_chr_map4.txt
    cut -f1-2 ${outfile%.*}_lifted_hg17_to_hg19.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map4.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg17_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg17_hg19.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map4.txt ${outfile%.*}_hg17_hg19.txt | cut -f1 > ${outfile%.*}_hg17_unmapped_b37.txt

    #liftover from hg17 to hg38
    $liftover ${outfile%.*}_hg17.txt $chainfile_hg1738 ${outfile%.*}_lifted_hg17_to_hg38.bed ${outfile%.*}_unlifted_hg17_to_hg38.bed
    rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map4.txt
    cut -f1-2 ${outfile%.*}_lifted_hg17_to_hg38.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map4.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg17_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg17_hg38.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map4.txt ${outfile%.*}_hg17_hg38.txt | cut -f1 > ${outfile%.*}_hg17_unmapped_b38.txt
fi

#chr:pos neither seen in both (check for hg16)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg17_unmapped_b38.txt ${outfile%.*}_hg17_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg16.txt
if [ -s ${outfile%.*}_hg16.txt ]; then
    #liftover from hg16 to hg19
    $liftover ${outfile%.*}_hg16.txt $chainfile_hg1619 ${outfile%.*}_lifted_hg16_to_hg19.bed ${outfile%.*}_unlifted_hg16_to_hg19.bed
    rm -f $outfolder/gwasc_${tag}_b37_dbsnp_chr_map5.txt
    cut -f1-2 ${outfile%.*}_lifted_hg16_to_hg19.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_${tag}_b37_dbsnp_chr_map5.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg16_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg16_hg19.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map5.txt ${outfile%.*}_hg16_hg19.txt | cut -f1 > ${outfile%.*}_hg16_unmapped_b37.txt

    #liftover from hg16 to hg38
    $liftover ${outfile%.*}_hg16.txt $chainfile_hg1638 ${outfile%.*}_lifted_hg16_to_hg38.bed ${outfile%.*}_unlifted_hg16_to_hg38.bed
    rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map5.txt
    cut -f1-2 ${outfile%.*}_lifted_hg16_to_hg38.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map5.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg16_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg16_hg38.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map5.txt ${outfile%.*}_hg16_hg38.txt | cut -f1 > ${outfile%.*}_hg16_unmapped_b38.txt
fi

#chr:pos neither seen in both (check for hg15)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg16_unmapped_b38.txt ${outfile%.*}_hg16_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg15.txt
if [ -s ${outfile%.*}_hg15.txt ]; then
    #liftover from hg15 to hg38
    $liftover ${outfile%.*}_hg15.txt $chainfile_hg1538 ${outfile%.*}_lifted_hg15_to_hg38.bed ${outfile%.*}_unlifted_hg15_to_hg38.bed
    rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map6.txt
    cut -f1-2 ${outfile%.*}_lifted_hg15_to_hg38.bed | sed 's/chr//' | while read line; do
        chr=`echo "$line" | cut -f1`
        snp=`echo "$line" | awk '{print $1":"$2}'`
        echo $chr $snp
        awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_${tag}_b38_dbsnp_chr_map6.txt
    done
    cut -f1-2,4 ${outfile%.*}_lifted_hg15_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' | sort -s -k2,2 > ${outfile%.*}_hg15_hg38.txt
    awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map6.txt ${outfile%.*}_hg15_hg38.txt | cut -f1 > ${outfile%.*}_hg15_unmapped_b38.txt
fi
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
awk -F '\t' '{print $1,$2,$2,$8,$9}' OFS='\t' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map0_2.txt
awk -F '\t' '{print $1,$2,$2,$8,$9}' OFS='\t' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0.txt > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map0_2.txt
for i in 18 17 16 15; do
    if [ $i -eq 18 ]; then
        j=3
    elif [ $i -eq 17 ]; then
        j=4
    elif [ $i -eq 16 ]; then
        j=5
    elif [ $i -eq 15 ]; then
        j=6
    fi

    awk -F '\t' '{print $1,$2,$8,$9}' OFS='\t' $outfolder/gwasc_${tag}_b38_dbsnp_chr_map${j}.txt | sort -s -k2,2 > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map${j}.tmp
    join -1 2 -2 2 -t '	' ${outfile%.*}_hg${i}_hg38.txt $outfolder/gwasc_${tag}_b38_dbsnp_chr_map${j}.tmp | awk -F '\t' '{print $3,$2,$1,$4,$5}' OFS='\t' > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map${j}_2.txt
    rm -f $outfolder/gwasc_${tag}_b38_dbsnp_chr_map${j}.tmp
    if [ $i -eq 15 ]; then
        continue
    fi

    awk -F '\t' '{print $1,$2,$8,$9}' OFS='\t' $outfolder/gwasc_${tag}_b37_dbsnp_chr_map${j}.txt | sort -s -k2,2 > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map${j}.tmp
    join -1 2 -2 2 -t '	' ${outfile%.*}_hg${i}_hg19.txt $outfolder/gwasc_${tag}_b37_dbsnp_chr_map${j}.tmp | awk -F '\t' '{print $3,$2,$1,$4,$5}' OFS='\t' > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map${j}_2.txt
    rm -f $outfolder/gwasc_${tag}_b37_dbsnp_chr_map${j}.tmp
done
fi

run_flag='Y'
if [ "$run_flag" = 'Y' ]; then
#build dbsnp rs map (output: rsid(1-based_coordinate) 0-based coordinate variant_type 1-based coordinate)
awk -F '\t' '{print $1,$9,$8,$2}' OFS='\t' $outfolder/gwasc_${tag}_b38_dbsnp_rs_map2.txt | sort | uniq | sort -s -k1,1 > $outfolder/gwasc_${tag}_b38_dbsnp_rs_map_2.txt
awk -F '\t' '{print $1,$9,$8,$2}' OFS='\t' $outfolder/gwasc_${tag}_b37_dbsnp_rs_map2.txt | sort | uniq | sort -s -k1,1 > $outfolder/gwasc_${tag}_b37_dbsnp_rs_map_2.txt
#build dbsnp chr map (output: rsid gwasc_0-based_coordinate 0-based coordinate variant_type 1-based coordinate)
cat $outfolder/gwasc_${tag}_b38_dbsnp_chr_map*_2.txt | sort | uniq | sort -s -k2,2 > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map_2.txt
cat $outfolder/gwasc_${tag}_b37_dbsnp_chr_map*_2.txt | sort | uniq | sort -s -k2,2 > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map_2.txt
#update gwasc ids according to maps
echo 'CHR	POS	CURRENT_ID	RSID	Hg19_ID	Hg38_ID	SNP	ID0	RISK_ALLELE	MERGED	ID	STUDY	TRAIT	FLAG	VAR_TYPE' > ${outfile%.*}_avail_2.txt
tail -n +2 ${outfile%.*}_avail.txt | sort -s -k3,3 > ${outfile%.*}_avail.tmp
#update gwasc rsid hg19_chr:pos in line with rs map
join -1 3 -2 1 -a 1 -t '	' ${outfile%.*}_avail.tmp $outfolder/gwasc_${tag}_b37_dbsnp_rs_map_2.txt | awk -F '\t' '{if(NF>11) print $2,$3,$1,$1,$12,$4,$5,$6,$7,$8,$9,$10,$11,$13; else print $2,$3,$1,$1,"-",$4,$5,$6,$7,$8,$9,$10,$11,"-"}' OFS='\t' | sort | uniq | sort -s -k3,3 > ${outfile%.*}_avail.tmp2
#update gwasc rsid hg38_chr:pos in line with rs map
join -1 3 -2 1 -a 1 -t '	' ${outfile%.*}_avail.tmp2 $outfolder/gwasc_${tag}_b38_dbsnp_rs_map_2.txt | awk -F '\t' '{if(NF>14) print $2,$3,$1,$4,$5,$15,$6,$7,$8,$9,$10,$11,$12,$13,$16; else print $2,$3,$1,$4,$5,"-",$6,$7,$8,$9,$10,$11,$12,$13,$14}' OFS='\t' | sort | uniq | sort -s -k3,3 > ${outfile%.*}_avail.tmp3
cut -f3,15 ${outfile%.*}_avail.tmp3 | grep rs | awk -F '\t' '$2=="-"' | cut -f1 | sort | uniq > ${outfile%.*}_rs_unfound.txt
#update gwasc rsid hg19_chr:pos in line with chr map
join -1 3 -2 2 -a 1 -t '	' ${outfile%.*}_avail.tmp3 $outfolder/gwasc_${tag}_b37_dbsnp_chr_map_2.txt | awk -F '\t' '{if(NF>15) print $2,$3,$1,$16,$17,$6,$7,$8,$9,$10,$11,$12,$13,$14,$18; else print $2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | sort | uniq | sort -s -k3,3 > ${outfile%.*}_avail.tmp4
#update gwasc rsid hg38_chr:pos in line with chr map
join -1 3 -2 2 -a 1 -t '	' ${outfile%.*}_avail.tmp4 $outfolder/gwasc_${tag}_b38_dbsnp_chr_map_2.txt | awk -F '\t' '{if(NF>15) print $2,$3,$1,$16,$5,$17,$7,$8,$9,$10,$11,$12,$13,$14,$18; else print $2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | sort | uniq | sort -s -k4,4 > ${outfile%.*}_avail.tmp5
awk -F '\t' '$4!~/rs/' ${outfile%.*}_avail.tmp5 | cut -f4 | sort | uniq > ${outfile%.*}_bp_unfound.txt
#refresh chr map rsid with both b37 and b38 coordinates
sort -s -k1,1 $outfolder/gwasc_${tag}_b37_dbsnp_chr_map_2.txt > $outfolder/gwasc_${tag}_b37_dbsnp_chr_map_2.tmp
sort -s -k1,1 $outfolder/gwasc_${tag}_b38_dbsnp_chr_map_2.txt > $outfolder/gwasc_${tag}_b38_dbsnp_chr_map_2.tmp
awk -F '\t' '$3~/:/' ${outfile%.*}_avail.tmp5 > ${outfile%.*}_avail.tmp51
join -1 4 -2 1 -a 1 -t '	' ${outfile%.*}_avail.tmp51 $outfolder/gwasc_${tag}_b37_dbsnp_chr_map_2.tmp | awk -F '\t' '{if(NF>15) print $2,$3,$4,$1,$17,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15; else print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | sort | uniq | sort -s -k4,4 > ${outfile%.*}_avail.tmp6
join -1 4 -2 1 -a 1 -t '	' ${outfile%.*}_avail.tmp6 $outfolder/gwasc_${tag}_b38_dbsnp_chr_map_2.tmp | awk -F '\t' '{if(NF>15) print $2,$3,$4,$1,$5,$17,$7,$8,$9,$10,$11,$12,$13,$14,$15; else print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' OFS='\t' | sort | uniq | sort -n -k1,1 -k2,2 > ${outfile%.*}_avail.tmp7
awk -F '\t' '$3!~/:/' ${outfile%.*}_avail.tmp5 >> ${outfile%.*}_avail.tmp7
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
echo "**********************************B38|B37**********************************"
echo "pos sum:"`cat ${outfile%.*}_bp.txt | wc -l`
n_summ=`cat $outfolder/gwasc_B_b3738_dbsnp_summ.txt | wc -l`
n_comm=`cat $outfolder/gwasc_B_b3738_dbsnp_comm.txt | wc -l`
n_miss=`cat $outfolder/gwasc_B_b3738_dbsnp_miss.txt | wc -l`
echo "map sum:"$n_summ" [common:"$n_comm"], miss:"$n_miss" (total:"$((n_summ+n_miss))")"
n38_1=`cut -f2 $outfolder/gwasc_B_b38_dbsnp_chr_map1.txt | sort | uniq | wc -l`
n38_2=`cut -f2 $outfolder/gwasc_B_b38_dbsnp_chr_map2.txt | sort | uniq | wc -l`
n37_2=`cut -f2 $outfolder/gwasc_B_b37_dbsnp_chr_map2.txt | sort | uniq | wc -l`
n38_0=`cut -f2 $outfolder/gwasc_B_b38_dbsnp_chr_map0.txt | sort | uniq | wc -l`
n37_0=`cut -f2 $outfolder/gwasc_B_b37_dbsnp_chr_map0.txt | sort | uniq | wc -l`
echo "map sum:"$((n38_1+n38_2+n37_2))" [b38_common:"$n38_1", b38_specific:"$n38_2", b37_specific:"$n37_2"]"
echo "map b38:"$n38_0" map b37:"$n37_0"]"

wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map1.txt
echo $n38_1"(uniq)"
head -n3 $outfolder/gwasc_B_b38_dbsnp_chr_map1.txt

wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map2.txt
echo $n38_2"(uniq)"
head -n3 $outfolder/gwasc_B_b38_dbsnp_chr_map2.txt

wc -l $outfolder/gwasc_B_b37_dbsnp_chr_map2.txt
echo $n37_2"(uniq)"
head -n3 $outfolder/gwasc_B_b37_dbsnp_chr_map2.txt

wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map0.txt
echo $n38_0"(uniq)"
head -n3 $outfolder/gwasc_B_b38_dbsnp_chr_map0.txt
wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map0_2.txt
head -n3 $outfolder/gwasc_B_b38_dbsnp_chr_map0_2.txt

wc -l $outfolder/gwasc_B_b37_dbsnp_chr_map0.txt
echo $n37_0"(uniq)"
head -n3 $outfolder/gwasc_B_b37_dbsnp_chr_map0.txt
wc -l $outfolder/gwasc_B_b37_dbsnp_chr_map0_2.txt
head -n3 $outfolder/gwasc_B_b37_dbsnp_chr_map0_2.txt

wc -l $outfolder/gwasc_B_b3738_dbsnp_miss.txt
head -n3 $outfolder/gwasc_B_b3738_dbsnp_miss.txt
echo ""

for i in 18 17 16 15; do
    echo "********************************* Hg$i ************************************"
    wc -l ${outfile%.*}_hg$i.txt
    head -n3 ${outfile%.*}_hg$i.txt
    echo "..."
    tail -n1 ${outfile%.*}_hg$i.txt
    echo ""

    if [ $i -eq 18 ]; then
        j=3
    elif [ $i -eq 17 ]; then
        j=4
    elif [ $i -eq 16 ]; then
        j=5
    elif [ $i -eq 15 ]; then
        j=6
    fi

    if [ $i -ne 15 ]; then
        echo "#----------------------liftover from Hg$i to Hg19-------------------------#"
        wc -l ${outfile%.*}_unlifted_hg${i}_to_hg19.bed | awk '{print $1/2, $2}'
        grep -v Del ${outfile%.*}_unlifted_hg${i}_to_hg19.bed | head -n2
        wc -l ${outfile%.*}_lifted_hg${i}_to_hg19.bed
        head -n10 ${outfile%.*}_lifted_hg${i}_to_hg19.bed
        wc -l ${outfile%.*}_hg${i}_hg19.txt
        head -n10 ${outfile%.*}_hg${i}_hg19.txt
        wc -l $outfolder/gwasc_B_b37_dbsnp_chr_map${j}.txt
        echo `cut -f2 $outfolder/gwasc_B_b37_dbsnp_chr_map${j}.txt | sort | uniq | wc -l`"(uniq bp)"
        head -n10 $outfolder/gwasc_B_b37_dbsnp_chr_map${j}.txt
        wc -l $outfolder/gwasc_B_b37_dbsnp_chr_map${j}_2.txt
        head -n10 $outfolder/gwasc_B_b37_dbsnp_chr_map${j}_2.txt
        wc -l ${outfile%.*}_hg${i}_unmapped_b37.txt
        n_hg19=`cat ${outfile%.*}_lifted_hg${i}_to_hg19.bed | wc -l`
        n_map=`cut -f2 $outfolder/gwasc_B_b37_dbsnp_chr_map${j}.txt | sort | uniq | wc -l`
        echo $((n_hg19-n_map))
        head -n10 ${outfile%.*}_hg${i}_unmapped_b37.txt
        echo ""
    fi

    echo "#----------------------liftover from Hg$i to Hg38-------------------------#"
    wc -l ${outfile%.*}_unlifted_hg${i}_to_hg38.bed | awk '{print $1/2, $2}'
    grep -v Del ${outfile%.*}_unlifted_hg${i}_to_hg38.bed | head -n2
    wc -l ${outfile%.*}_lifted_hg${i}_to_hg38.bed
    head -n10 ${outfile%.*}_lifted_hg${i}_to_hg38.bed
    wc -l ${outfile%.*}_hg${i}_hg38.txt
    head -n10 ${outfile%.*}_hg${i}_hg38.txt
    wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map${j}.txt
    echo `cut -f2 $outfolder/gwasc_B_b38_dbsnp_chr_map${j}.txt | sort | uniq | wc -l`"(uniq bp)"
    head -n10 $outfolder/gwasc_B_b38_dbsnp_chr_map${j}.txt
    wc -l $outfolder/gwasc_B_b38_dbsnp_chr_map${j}_2.txt
    head -n10 $outfolder/gwasc_B_b38_dbsnp_chr_map${j}_2.txt
    wc -l ${outfile%.*}_hg${i}_unmapped_b38.txt
    n_hg38=`cat ${outfile%.*}_lifted_hg${i}_to_hg38.bed | wc -l`
    n_map=`cut -f2 $outfolder/gwasc_B_b38_dbsnp_chr_map${j}.txt | sort | uniq | wc -l`
    echo $((n_hg38-n_map))
    head -n10 ${outfile%.*}_hg${i}_unmapped_b38.txt
    echo ""
done
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build b38 and b37 map for snps with only chr:pos by online search
echo $outfolder
rm -f $outfolder/dbsnp_b38_base_1.txt
rm -f $outfolder/dbsnp_b37_base_1.txt
#for id in `cat ${outfile%.*}_bp.txt`; do
for id in `cat ${outfile%.*}_hg15_unmapped_b38.txt`; do
    #id=4:74265673 #id=8:95855720 #id=14:69912282 #id=17:7536691 #id=X:65638484 id=7:2060397 id=9:31295759
    echo $id
    chr=${id%:*}
    pos=${id#*:}
    mkdir -p $dbsnp_map_folder/web
    if [ ! -f $dbsnp_map_folder/web/term=$id ]; then
        wget https://www.ncbi.nlm.nih.gov/snp/?term=${chr}%3A${pos} -O $dbsnp_map_folder/web/term=$id
    fi
    if [ -s $dbsnp_map_folder/web/term=$id ] && [ `grep Chromosome: $dbsnp_map_folder/web/term=$id | grep -v 'no mapping' | wc -l` -gt 0 ]; then
        awk '/has merged into/{f=0} f; /begin Results/{f=1}' $dbsnp_map_folder/web/term=$id > $outfolder/term=$id.tmp
        p38=`grep Chromosome: $outfolder/term=$id.tmp | grep $pos | head -n1 | sed 's/.*background-color:">//' | sed 's/<\/span>//' | sed 's/.*://'`
        p37=`grep GRCh38 $outfolder/term=$id.tmp | grep $pos |  head -n1 | sed 's/.*background-color:">//' | sed 's/<\/span>//' | sed 's/.*://'`
        if [ "$p38" != "" ]; then
            snp=`grep -B 2 Chromosome: $outfolder/term=$id.tmp | grep -B 2 \>$pos | grep 'Homo sapiens' | grep -v 'has merged into' | sed 's/.*snp_info="//' | cut -d: -f1`
            p37=`grep -A 1 Chromosome: $outfolder/term=$id.tmp | grep -A 1 \>$pos | grep GRCh | grep -v \>$pos\< | sed 's/.*background-color:">//' | sed 's/<\/span>//' | sed 's/.*://' | uniq`
        
        else
            snp=`grep -B 3 GRCh38 $outfolder/term=$id.tmp | grep -B 3 \>$pos | grep 'Homo sapiens' | grep -v 'has merged into' | sed 's/.*snp_info="//' | cut -d: -f1`
            p38=`grep -B 1 GRCh38 $outfolder/term=$id.tmp | grep -B 1 \>$pos | grep Chromosome: | grep -v \>$pos\< | sed 's/.*background-color:">//' | sed 's/<\/span>//' | sed 's/.*://' | uniq`
        fi
        if [ `echo $snp | sed 's/ /\n/g' | wc -l` -gt 1 ]; then
            echo "SKIP: more than one variant found for id=$id:" $snp
            rm -f $outfolder/term=$id.tmp
            continue
        fi
        if [ ! -f $dbsnp_map_folder/web/$snp ]; then
            wget https://www.ncbi.nlm.nih.gov/snp/$snp -O $dbsnp_map_folder/web/$snp
        fi
        if [ -s $dbsnp_map_folder/web/$snp ]; then
            a38=`cat $dbsnp_map_folder/web/$snp | grep -A50 'Genomic Placements' | grep GRCh38 | cut -d'>' -f2 | cut -d' ' -f1 | sort | tail -n1`
            a37=`cat $dbsnp_map_folder/web/$snp | grep -A50 'Genomic Placements' | grep GRCh37 | cut -d'>' -f2 | cut -d' ' -f1 | sort | tail -n1`
            s38=`cat $dbsnp_map_folder/web/$snp | grep -A50 'Genomic Placements' | grep -A1 GRCh38 | grep -v GRCh | grep td | cut -d'>' -f2 | cut -d: -f1 | sort | tail -n1`
            s37=`cat $dbsnp_map_folder/web/$snp | grep -A50 'Genomic Placements' | grep -A1 GRCh37 | grep -v GRCh | grep td | cut -d'>' -f2 | cut -d: -f1 | sort | tail -n1`
            echo "$snp	$chr:$p38	$chr	$p38	$a38	$s38" >> $outfolder/dbsnp_b38_base_1.txt
            echo "$snp	$chr:$p37	$chr	$p37	$a37	$s37" >> $outfolder/dbsnp_b37_base_1.txt
        fi
        rm -f $outfolder/term=$id.tmp
    fi
done
fi
