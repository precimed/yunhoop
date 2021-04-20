#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to ensemble NHGRI-EBI GWAS Catalog to make a clean file
# for novelty checking of loci WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 3 ]; then
  echo "Usage: sh compile_gwasc.sh gwasc build_id outfile"
  echo "Arguments: gwasc - gwascatalog"
  echo "           dbsnp_map_folder - top folder to hold built dbsnp maps within json and vcf subfolders"
  echo "           outfile - output file"
  echo "Example: sh compile_gwasc.sh gwas_catalog_v1.0-associations_e100_r2021-04-12.tsv $lf/ncbi/b153 gwas_catalog_v1.0-associations_e100_r2021-04-12.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
dbsnp_map_folder=$2
outfile=$3

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

#process records with only single or multiple ids (chr:pos or rs)
echo "process records with only single or multiple ids"
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/' > $gwasc.tmp1
rm -f $gwasc.tmp10
cut -f3- $gwasc.tmp1 | while read line; do
    rest=`echo "$line" | cut -f2-5`
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
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/ && $3!~/x/' $gwasc.tmp | tail -n +2 | cut -f1-5 > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/del-//gi'| sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f3 $gwasc.tmp2 > $gwasc.tmp23
cut -f4-7 $gwasc.tmp2 > $gwasc.tmp24
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 $gwasc.tmp24 > $gwasc.tmp2

#process records associated with multiple snps with chr:pos and ids
echo "process records with multiple snps with chr:pos and ids"
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/ || $3~/x/)' $gwasc.tmp > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp31
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp32
    echo "$line" | cut -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
    rest=`echo "$line" | cut -f4-7`
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
cut -f3 $outfile | grep rs | sort | uniq > ${outfile%.*}_rs.txt
fi

run_flag='Y'
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
if [ ! -f $dbsnp_map_folder/vcf/dbsnp_b38_chr1.txt]; then
fi

#build gwasc snp merge map (based on gwasc rs list and dbsnp_merge_map.txt)
sh $(dirname $0)/build_snp_merge_map.sh $dbsnp_map_folder/json/dbsnp_merge_map.txt ${outfile%.*}_rs.txt $(dirname $outfile)/gwasc_snp_merge_map.txt

#update gwasc rs list with merged-into snps (based on gwasc_snp_merge_map.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $(dirname $outfile)/gwasc_snp_merge_map.txt ${outfile%.*}_rs.txt > ${outfile%.*}_rs_merged.txt
cut -f2 $(dirname $outfile)/gwasc_snp_merge_map.txt | sort | uniq >> ${outfile%.*}_rs_merged.txt

#divide gwasc rs list into available and unavailable rs list (based on updated rs list and dbsnp_unavailable.txt)
awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $dbsnp_map_folder/json/dbsnp_unavailable.txt ${outfile%.*}_rs_merged.txt > ${outfile%.*}_rs_avail.txt
diff ${outfile%.*}_rs_merged.txt ${outfile%.*}_rs_avail.txt | grep \< | awk '{print $2}' > ${outfile%.*}_rs_unavail.txt

#build gwasc snp map (based on available rs list and )
sh $(dirname $0)/build_snp_map.sh $dbsnp_map_folder ${outfile%.*}_rs_avail.txt $(dirname $outfile)

awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/gwasc_snp_merge_map.txt $outfile > ${outfile%.*}_merged.txt
awk -F '\t' 'NR==FNR {D[$1]++; next} ($3 in D)' $(dirname $outfile)/gwasc_snp_merge_map.txt $outfile | sort -s -k3,3 > ${outfile%.*}.tmp
sort -s -k1,1 $(dirname $outfile)/gwasc_snp_merge_map.txt > $(dirname $outfile)/gwasc_snp_merge_map.tmp
join -1 3 -2 1 -t '	' ${outfile%.*}.tmp $(dirname $outfile)/gwasc_snp_merge_map.tmp | awk -F '\t' '{if(NF>9) print $2,$3,$10,$1,$4,$5,$6,$7,$8,$9; else print $2,$3,$1,$1,$4,$5,$6,$7,$8,$9}' OFS='\t' > ${outfile%.*}.tmp2
fi
