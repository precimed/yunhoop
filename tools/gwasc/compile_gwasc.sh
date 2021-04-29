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
  echo "           reffile - reference file"
  echo "           outfile - output file"
  echo "Example: sh compile_gwasc.sh gwas_catalog_v1.0-associations_e100_r2021-04-12.tsv $lf/ncbi/b153 $data/ref/9545380.ref gwas_catalog_v1.0-associations_e100_r2021-04-12.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
dbsnp_map_folder=$2
reffile=$3
outfile=$4
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

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
cut -f7-8,12-13,22-24 $gwasc | awk -F '\t' '{print $3,$4,$5,$6,$7,$1,$2}' OFS='\t' | tr -dc '\0-\177' | sed 's/&beta;/Beta /g' > $gwasc.txt

#filter snps with irregular ids
echo "filter irregular snps with empty chr"
cut -f1-3 $gwasc.txt | awk -F '\t' '$1==""' | grep -v ":" | grep -i -v rs | grep -i -v chr | cut -f3 > $gwasc.tmp
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $gwasc.tmp $gwasc.txt > $gwasc.tmp2
mv $gwasc.tmp2 $gwasc.txt

#clean snp ids
echo "clean snps"
head -n1 $gwasc.txt > $gwasc.tmp
tail -n +2 $gwasc.txt | cut -f1-4 | sed 's/\./:/' | sed 's/\s+//g' | sed 's/ ; /;/g' | sed 's/; /;/g' | sed 's/ ;/;/g' | sed 's/ , /;/g' | sed 's/, /;/g' | sed 's/ ,/;/g' | sed 's/,/;/g' | sed 's/ \/ /,/g' | sed 's/ \//,/g' | sed 's/\/ /,/g' | sed 's/\//,/g' | sed 's/ x /x/gi' > $gwasc.tmp1
tail -n +2 $gwasc.txt | cut -f5 | cut -d'_' -f1 | cut -d'-' -f1 | sed 's/[a-zA-Z]$//' | sed 's/^/rs/' | sed 's/^rs$/-/' > $gwasc.tmp2
tail -n +2 $gwasc.txt | cut -f6-7 > $gwasc.tmp3
for((i=0; i<=22; i++)); do
    sed "s/chr${i}: /chr${i}:/gi" $gwasc.tmp1 > $gwasc.tmp10
    sed "s/chr${i}_/chr${i}:/gi" $gwasc.tmp10 > $gwasc.tmp1
done
sed 's/chr://gi' $gwasc.tmp1 | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/		.*psy_/		/gi' | sed 's/hg18_/hg18/gi' > $gwasc.tmp10
paste -d'	' $gwasc.tmp10 $gwasc.tmp2 $gwasc.tmp3 >> $gwasc.tmp
rm -f $gwasc.tmp1* $gwasc.tmp2 $gwasc.tmp3

if [ `awk -F '\t' '($3~/x/ || $3~/;/)&&$5!="-" {print $3,$5}' OFS='\t' $gwasc.tmp | wc -l` -gt 0 ]; then
    echo "multi ids with single current id"
    awk -F '\t' '($3~/x/ || $3~/;/)&&$5!="-" {print $3,$5}' OFS='\t' $gwasc.tmp
    exit 1
fi

#process records with only single or multiple ids (chr:pos or rs)
echo "process records with only single or multiple ids"
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/' > $gwasc.tmp1
rm -f $gwasc.tmp10
cut -f3- $gwasc.tmp1 | while read line; do
    rest=`echo "$line" | cut -f2-`
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
done
mv $gwasc.tmp10 $gwasc.tmp1

#process records as single snps with chr:pos and id
echo "process records as single snp with chr:pos and id (chr:pos or rs)"
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/ && $3!~/x/' $gwasc.tmp | tail -n +2 > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/del-//gi'| sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f3 $gwasc.tmp2 > $gwasc.tmp23
cut -f4- $gwasc.tmp2 > $gwasc.tmp24
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 $gwasc.tmp24 > $gwasc.tmp2

#process records associated with multiple snps with chr:pos and ids
echo "process records with multiple snps with chr:pos and ids"
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/ || $3~/x/)' $gwasc.tmp > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp31
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp32
    echo "$line" | cut -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
    rest=`echo "$line" | cut -f4-`
    n1=`cat $gwasc.tmp31 | wc -l`
    n2=`cat $gwasc.tmp32 | wc -l`
    n3=`cat $gwasc.tmp33 | wc -l`
    if [ $n1 -ne $n2 ]; then
        echo $line
        echo "NOTE: chr number ($n1) does not match pos number ($n2)"
        exit 1
    fi
    if [ $n1 -gt $n3 ]; then
        echo $line
        echo "NOTE: chr number ($n1) is greater than snp number ($n3)"
        exit 1
    fi
    # a) number of chr:pos matches number of ids
    if [ $n3 -eq $n1 ] && [ `echo "$line" | cut -f3 | grep 'x' | wc -l` -gt 0 ]; then
        rm -f $gwasc.tmp34
        for ((i=1; i<=$n3; i++)); do
            echo "$rest" >> $gwasc.tmp34
        done
    # b) single chr:pos with multiple merging rs
    elif [ $n1 -eq 1 ] && [ $n3 -gt $n1 ] && [ `echo "$line" | cut -f3 | grep ',' | wc -l` -gt 0 ]; then
        for((i=2; i<=$n3; i++)); do
            head -n1 $gwasc.tmp31 >> $gwasc.tmp31
            head -n1 $gwasc.tmp32 >> $gwasc.tmp32
        done
        rm -f $gwasc.tmp34
        for ((i=1; i<=$n3; i++)); do
            echo "$rest" >> $gwasc.tmp34
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
        rm -f $gwasc.tmp34
        for((i=1; i<=$n3; i++)); do
            echo '-' >> $gwasc.tmp31
            echo '-' >> $gwasc.tmp32
            echo "$rest" >> $gwasc.tmp34
        done
    fi
    sed 's/del-//gi' $gwasc.tmp33 | sed 's/[a-zA-Z]$//' > $gwasc.tmp35
    paste -d'	' $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp35 $gwasc.tmp33 $gwasc.tmp34 >> $gwasc.tmp30
    rm -f $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp35 $gwasc.tmp33 $gwasc.tmp34
done
mv $gwasc.tmp30 $gwasc.tmp3

cat $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3 | sed 's/^23	/X	/' | sed 's/	23:/	X:/' | sed 's/^24	/Y	/' | sed 's/	24:/	Y:/' | sort | uniq > $outfile

#update chr:pos of records with only rs in light of available snps
echo "update chr:pos of records with only rs in light of available snps"
tail -n +2 $outfile | grep -v ^- > $outfile.tmp
cut -f1-3 $outfile.tmp | sort -s -k3,3 | uniq > $outfile.tmp1
grep ^- $outfile | sort -s -k3,3 > $outfile.tmp2
join -1 3 -2 3 -a 1 -t '	' $outfile.tmp2 $outfile.tmp1 | awk -F '\t' '{if (NF>8) print $9,$10,$1,$4,$5,$6,$7,$8; else print $2,$3,$1,$4,$5,$6,$7,$8}' OFS='\t' >> $outfile.tmp

#create initial set (OUTPUT: CHR POS SNP ID STUDY TRAIT)
echo "compile initial set"
head -n1 $outfile > $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1~/^[0-9]+$/' | sort -n -k1,1 -k2,2 >> $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1!~/^[0-9]+$/' | sort -k1,1 >> $outfile.tmp2
grep ^- $outfile.tmp >> $outfile.tmp2
echo 'CHR	POS	SNP	ID	MERGED	CURRENT_ID	STUDY	TRAIT	FLAG' > $outfile
#flag: 0-snp with only rs; 1-snp with chr:pos; 2-snp with both rs and chr:pos
awk -F '\t' '{if($1=="-") print $1,$2,$3,$4,$5,$6,$7,$8,"0"; else if($1!="-"&&$3~/:/) print $1,$2,$3,$4,$5,$6,$7,$8,"1"; else if($1!="-"&&$3~/rs/) print $1,$2,$3,$4,$5,$6,$7,$8,"2"; else print $1,$2,$3,$4,$5,$6,$7,$8,"x"}' OFS='\t' $outfile.tmp2 >> $outfile

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
rm -f $gwasc.tmp*
rm -f $outfile.tmp*

awk -F '\t' '{if($6=="-") print $1,$2,$3,$4,$5,$3,$7,$8,$9; else print $0}' OFS='\t' $outfile > ${outfile%.*}_cur.csv
cut -f6 ${outfile%.*}_cur.csv | grep rs | sort | uniq > ${outfile%.*}_rs.txt
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#chunksize=15000
#awk -F '\t' '$7=="0" || $7=="2"' $outfile | cut -f3 | sort | uniq > ${outfile%.*}_A.txt
#split -l $chunksize --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".txt" ${outfile%.*}_A.txt ${outfile%.*}_A_

#for i in ${outfile%.*}_A_*.txt; do
#    suffix=${i%.*}
#    suffix=${suffix##*_}
#    sh $(dirname $0)/search_dbsnp.sh $i ../dbsnp . $suffix &> ${i%.*}.log &
#done

#build dbSNP misc maps
if [ ! -f $dbsnp_map_folder/json/dbsnp_merge_map.txt ]; then
    sh $(dirname $0)/../dbsnp/build_dbsnp_maps.sh $dbsnp_map_folder N Y
    cat $dbsnp_map_folder/json/dbsnp_unsupported.txt $dbsnp_map_folder/json/dbsnp_withdrawn.txt $dbsnp_map_folder/json/dbsnp_nosnppos.txt > $dbsnp_map_folder/json/dbsnp_unavailable.txt
fi

#build dbSNP chr maps
if [ ! -f $dbsnp_map_folder/vcf/dbsnp_b38_chr1.txt ]; then
    sh $(dirname $0)/../dbsnp/build_dbsnp_maps.sh $dbsnp_map_folder Y N
fi

#build ref snp merge map (based on ref rs list [9545380.ref] and dbsnp_merge_map.txt)
if [ ! -f ${reffile%.*}_snp_merge_map.txt ]; then
    tail -n +2 $reffile | cut -f2 > ${reffile%.*}_rs.txt
    sh $(dirname $0)/build_snp_merge_map.sh $dbsnp_map_folder/json/dbsnp_merge_map.txt ${reffile%.*}_rs.txt ${reffile%.*}_snp_merge_map.txt
    sort -s -k2,2 $reffile > ${reffile%.*}.sorted
    sort -s -k1,1 ${reffile%.*}_snp_merge_map.txt > ${reffile%.*}_snp_merge_map.sorted
    echo 'CHR	CURRENT_ID	SNP	BP	A1	A2	MERGED' > ${reffile%.*}.txt
    join -1 2 -2 1 -a 1 ${reffile%.*}.sorted ${reffile%.*}_snp_merge_map.sorted -t '	' | grep -v SNP | awk -F '\t' '{if(NF>6) print $2,$7,$1,$4,$5,$6,$8; else print $2,$1,$1,$4,$5,$6,"0"}' OFS='\t' | sort -n -k1,1 -k4,4 >> ${reffile%.*}.txt
    rm -f ${reffile%.*}.sorted ${reffile%.*}_snp_merge_map.sorted
fi

#get updated dbsnp b37/b38 chr:pos after merging
echo 'CURRENT_RS	CHR	SNP	BP	A1	A2	MERGED	POS	REF	ALT	ID	CURRENT_ID' > ${reffile%.*}_b37.txt
echo 'CURRENT_RS	CHR	SNP	BP	A1	A2	MERGED	POS	REF	ALT	ID	CURRENT_ID' > ${reffile%.*}_b38.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    echo $i
    grep "^$i	" ${reffile%.*}.txt | sort -s -k2,2 > ${reffile%.*}_chr$i.tmp
    awk -F '\t' 'NR==FNR {D[$2]++; next} ($1 in D)' ${reffile%.*}_chr$i.tmp $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt | sort -s -k1,1 > $(dirname $reffile)/dbsnp_b37_chr$i.tmp
    awk -F '\t' 'NR==FNR {D[$2]++; next} ($1 in D)' ${reffile%.*}_chr$i.tmp $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt | sort -s -k1,1 > $(dirname $reffile)/dbsnp_b38_chr$i.tmp
    join -1 2 -2 1 ${reffile%.*}_chr$i.tmp $(dirname $reffile)/dbsnp_b37_chr$i.tmp -t '	' | cut -f1-7,10-12 | awk -F '\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2":"$4,$2":"$8}' OFS='\t' | sort -n -k2,2 -k4,4 >> ${reffile%.*}_b37.txt
    join -1 2 -2 1 ${reffile%.*}_chr$i.tmp $(dirname $reffile)/dbsnp_b38_chr$i.tmp -t '	' | cut -f1-7,10-12 | awk -F '\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2":"$4,$2":"$8}' OFS='\t' | sort -n -k2,2 -k4,4 >> ${reffile%.*}_b38.txt
done
rm -f ${reffile%.*}_chr*.tmp $(dirname $reffile)/dbsnp_b*_chr*.tmp

#build gwasc snp merge map (based on gwasc rs list and dbsnp_merge_map.txt)
sh $(dirname $0)/build_snp_merge_map.sh $dbsnp_map_folder/json/dbsnp_merge_map.txt ${outfile%.*}_rs.txt $outfolder/gwasc_snp_merge_map.txt

#update gwasc rs list with merged-into snps (based on gwasc_snp_merge_map.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_rs.txt > ${outfile%.*}_rs_merged.txt
cut -f2 $outfolder/gwasc_snp_merge_map.txt | sort | uniq >> ${outfile%.*}_rs_merged.txt

#update gwasc rs list with available snps (based on dbsnp_unavailable.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_rs_merged.txt | sort | uniq > ${outfile%.*}_rs_avail.txt

#build gwasc snp chr map (based on available rs list and dbsnp chr map)
sh $(dirname $0)/build_snp_map.sh $dbsnp_map_folder ${outfile%.*}_rs_avail.txt $outfolder

#update in light of gwasc snp merge map
awk -F '\t' 'NR==FNR {D[$1]++; next} !($6 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_cur.csv | awk -F '\t' '{print $1,$2,$6,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' > ${outfile%.*}_merged.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($6 in D)' $outfolder/gwasc_snp_merge_map.txt ${outfile%.*}_cur.csv | sort -s -k6,6 > ${outfile%.*}.tmp
sort -s -k1,1 $outfolder/gwasc_snp_merge_map.txt > $outfolder/gwasc_snp_merge_map.tmp
join -1 6 -2 1 -a 1 -t '	' ${outfile%.*}.tmp $outfolder/gwasc_snp_merge_map.tmp | awk -F '\t' '{if(NF>9) print $2,$3,$10,$4,$5,$6,$1,$7,$8,$9; else print $2,$3,$1,$4,$5,$6,$1,$7,$8,$9}' OFS='\t' >> ${outfile%.*}_merged.txt
rm -f ${outfile%.*}.tmp $outfolder/gwasc_snp_merge_map.tmp

#update in light of dbsnp_unavailable.txt
head -n1 ${outfile%.*}_merged.txt > ${outfile%.*}_avail.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_merged.txt | tail -n +2 | sort | uniq >> ${outfile%.*}_avail.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($3 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_merged.txt | cut -f3 | sort | uniq > ${outfile%.*}_unavail.txt

if [ `cut -f3,7 ${outfile%.*}_avail.txt | sort | uniq | cut -f2 | sort | uniq -c | awk '$1>1' | wc -l` -gt 0 ]; then
    echo "snps merged more than once:"
    cut -f3,7 ${outfile%.*}_avail.txt | sort | uniq | cut -f2 | sort | uniq -c | awk '$1>1'
    exit 1
fi

#snps with only chr:pos
cut -f3 ${outfile%.*}_avail.txt | grep : | grep -v hg |  sort -s | uniq > ${outfile%.*}_bp.txt
fi

#build b38 and b37 map for snps with rs
rm -f $outfolder/gwasc_A_b38_dbsnp_rs_map.txt
rm -f $outfolder/gwasc_A_b37_dbsnp_rs_map.txt
rm -f $outfolder/gwasc_B_b38_dbsnp_rs_map.txt
rm -f $outfolder/gwasc_B_b37_dbsnp_rs_map.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    echo $i
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b38.txt >> $outfolder/gwasc_A_b38_dbsnp_rs_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b37.txt >> $outfolder/gwasc_A_b37_dbsnp_rs_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt >> $outfolder/gwasc_B_b38_dbsnp_rs_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_rs_avail.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt >> $outfolder/gwasc_B_b37_dbsnp_rs_map.txt
done

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
#build b38 and b37 map for snps with only chr:pos
rm -f $outfolder/gwasc_A_b38_dbsnp_snv_map.txt
rm -f $outfolder/gwasc_A_b37_dbsnp_snv_map.txt
rm -f $outfolder/gwasc_B_b38_dbsnp_snv_map.txt
rm -f $outfolder/gwasc_B_b37_dbsnp_snv_map.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    echo $i
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b38.txt | awk -F '\t' '$10=="SNV"' >> $outfolder/gwasc_A_b38_dbsnp_snv_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/json/dbsnp_chr${i}_b37.txt | awk -F '\t' '$10=="SNV"' >> $outfolder/gwasc_A_b37_dbsnp_snv_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b38_chr$i.txt | awk -F '\t' '$8=="SNV"'>> $outfolder/gwasc_B_b38_dbsnp_snv_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' ${outfile%.*}_bp.txt $dbsnp_map_folder/vcf/dbsnp_b37_chr$i.txt | awk -F '\t' '$8=="SNV"' >> $outfolder/gwasc_B_b37_dbsnp_snv_map.txt
done

#chr:pos shared between b38 and b37 are taken as b38 (since it currently mapped to genome assembly GRCh38)
awk -F '\t' 'NR==FNR {D[$2]++; next} ($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map.txt $outfolder/gwasc_B_b37_dbsnp_chr_map.txt | cut -f2| sort | uniq > $outfolder/gwasc_B_b3738_dbsnp_comm.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} ($2 in D)' $outfolder/gwasc_B_b37_dbsnp_chr_map.txt $outfolder/gwasc_B_b38_dbsnp_chr_map.txt > $outfolder/gwasc_B_b38_dbsnp_chr_map1.txt

#chr:pos only seen with b38 taken as b38
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b37_dbsnp_chr_map.txt $outfolder/gwasc_B_b38_dbsnp_chr_map.txt > $outfolder/gwasc_B_b38_dbsnp_chr_map2.txt

#chr:pos only seen with b37 taken as b37
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map.txt $outfolder/gwasc_B_b37_dbsnp_chr_map.txt > $outfolder/gwasc_B_b37_dbsnp_chr_map2.txt

#chr:pos not seen in both (check for hg18)
cat $outfolder/gwasc_B_b37_dbsnp_chr_map.txt $outfolder/gwasc_B_b38_dbsnp_chr_map.txt | cut -f2 | sort | uniq > $outfolder/gwasc_B_b3738_dbsnp_summ.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $outfolder/gwasc_B_b3738_dbsnp_summ.txt ${outfile%.*}_bp.txt | sort | uniq > $outfolder/gwasc_B_b3738_dbsnp_miss.txt
cat $outfolder/gwasc_B_b3738_dbsnp_miss.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg18.txt
cut -f3 ${outfile%.*}_avail.txt | grep hg18 | sed 's/hg18//' | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' >> ${outfile%.*}_hg18.txt
#liftover from hg18 to hg19
$liftover ${outfile%.*}_hg18.txt $chainfile_hg1819 ${outfile%.*}_lifted_hg18_to_hg19.bed ${outfile%.*}_unlifted_hg18_to_hg19.bed
rm -f $outfolder/gwasc_B_b37_dbsnp_chr_map3.txt
cut -f1-2 ${outfile%.*}_lifted_hg18_to_hg19.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_B_b37_dbsnp_chr_map3.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg18_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg18_hg19.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b37_dbsnp_chr_map3.txt ${outfile%.*}_hg18_hg19.txt | cut -f1 > ${outfile%.*}_hg18_unmapped_b37.txt

#liftover from hg18 to hg38
$liftover ${outfile%.*}_hg18.txt $chainfile_hg1838 ${outfile%.*}_lifted_hg18_to_hg38.bed ${outfile%.*}_unlifted_hg18_to_hg38.bed
rm -f $outfolder/gwasc_B_b38_dbsnp_chr_map3.txt
cut -f1-2 ${outfile%.*}_lifted_hg18_to_hg38.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_B_b38_dbsnp_chr_map3.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg18_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg18_hg38.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map3.txt ${outfile%.*}_hg18_hg38.txt | cut -f1 > ${outfile%.*}_hg18_unmapped_b38.txt

#chr:pos neither seen in both (check for hg17)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg18_unmapped_b38.txt ${outfile%.*}_hg18_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg17.txt
#liftover from hg17 to hg19
$liftover ${outfile%.*}_hg17.txt $chainfile_hg1719 ${outfile%.*}_lifted_hg17_to_hg19.bed ${outfile%.*}_unlifted_hg17_to_hg19.bed
rm -f $outfolder/gwasc_B_b37_dbsnp_chr_map4.txt
cut -f1-2 ${outfile%.*}_lifted_hg17_to_hg19.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_B_b37_dbsnp_chr_map4.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg17_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg17_hg19.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b37_dbsnp_chr_map4.txt ${outfile%.*}_hg17_hg19.txt | cut -f1 > ${outfile%.*}_hg17_unmapped_b37.txt
#liftover from hg17 to hg38
$liftover ${outfile%.*}_hg17.txt $chainfile_hg1738 ${outfile%.*}_lifted_hg17_to_hg38.bed ${outfile%.*}_unlifted_hg17_to_hg38.bed
rm -f $outfolder/gwasc_B_b38_dbsnp_chr_map4.txt
cut -f1-2 ${outfile%.*}_lifted_hg17_to_hg38.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_B_b38_dbsnp_chr_map4.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg17_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg17_hg38.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map4.txt ${outfile%.*}_hg17_hg38.txt | cut -f1 > ${outfile%.*}_hg17_unmapped_b38.txt

#chr:pos neither seen in both (check for hg16)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg17_unmapped_b38.txt ${outfile%.*}_hg17_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg16.txt
#liftover from hg16 to hg19
$liftover ${outfile%.*}_hg16.txt $chainfile_hg1619 ${outfile%.*}_lifted_hg16_to_hg19.bed ${outfile%.*}_unlifted_hg16_to_hg19.bed
rm -f $outfolder/gwasc_B_b37_dbsnp_chr_map5.txt
cut -f1-2 ${outfile%.*}_lifted_hg16_to_hg19.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b37_chr$chr.txt >> $outfolder/gwasc_B_b37_dbsnp_chr_map5.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg16_to_hg19.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg16_hg19.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b37_dbsnp_chr_map5.txt ${outfile%.*}_hg16_hg19.txt | cut -f1 > ${outfile%.*}_hg16_unmapped_b37.txt
#liftover from hg16 to hg38
$liftover ${outfile%.*}_hg16.txt $chainfile_hg1638 ${outfile%.*}_lifted_hg16_to_hg38.bed ${outfile%.*}_unlifted_hg16_to_hg38.bed
rm -f $outfolder/gwasc_B_b38_dbsnp_chr_map5.txt
cut -f1-2 ${outfile%.*}_lifted_hg16_to_hg38.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_B_b38_dbsnp_chr_map5.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg16_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg16_hg38.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map5.txt ${outfile%.*}_hg16_hg38.txt | cut -f1 > ${outfile%.*}_hg16_unmapped_b38.txt

#chr:pos neither seen in both (check for hg15)
awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' ${outfile%.*}_hg16_unmapped_b38.txt ${outfile%.*}_hg16_unmapped_b37.txt | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}_hg15.txt
#liftover from hg15 to hg38
$liftover ${outfile%.*}_hg15.txt $chainfile_hg1538 ${outfile%.*}_lifted_hg15_to_hg38.bed ${outfile%.*}_unlifted_hg15_to_hg38.bed
rm -f $outfolder/gwasc_B_b38_dbsnp_chr_map6.txt
cut -f1-2 ${outfile%.*}_lifted_hg15_to_hg38.bed | sed 's/chr//' | while read line; do
    chr=`echo "$line" | cut -f1`
    snp=`echo "$line" | awk '{print $1":"$2}'`
    echo $chr $snp
    awk -v snp=$snp -F '\t' '$2==snp' $dbsnp_map_folder/vcf/dbsnp_b38_chr$chr.txt >> $outfolder/gwasc_B_b38_dbsnp_chr_map6.txt
done
cut -f1-2,4 ${outfile%.*}_lifted_hg15_to_hg38.bed | sed 's/chr//' |  awk '{print $3,$1":"$2}' OFS='\t' > ${outfile%.*}_hg15_hg38.txt
awk -F '\t' 'NR==FNR {D[$2]++; next} !($2 in D)' $outfolder/gwasc_B_b38_dbsnp_chr_map6.txt ${outfile%.*}_hg15_hg38.txt | cut -f1 > ${outfile%.*}_hg15_unmapped_b38.txt
fi

run_flag='N'
if [ "$run_flag" = 'Y' ]; then
echo $outfolder
rm -f $outfolder/dbsnp_b38_base_1.txt
rm -f $outfolder/dbsnp_b37_base_1.txt
#for id in `cat $outfolder/gwasc_B_b3738_dbsnp_miss.txt`; do
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
