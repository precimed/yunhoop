#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to reform NHGRI-EBI GWAS Catalog to make a clean file
# for novelty checking of loci WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -ne 2 ]; then
  echo "Usage: sh reform_gwasc.sh gwasc outfile"
  echo "Arguments: gwasc - gwascatalog"
  echo "           outfile - output file"
  echo "Example: sh reform_gwasc22.sh gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv gwas_catalog_v1.0-associations_e100_r2021-02-25.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
outfile=$2

build37_snpchrpos_ref=$yc/software/liftover/ncbi/b132_SNPChrPosOnRef_37_1.bcp.gz
hg19_chrpossnp_ref=$yc/software/liftover/ucsc/hg19/hg19_snp151.txt.gz
hg38_chrpossnp_ref=$yc/software/liftover/ucsc/hg38/hg38_snp151.txt.gz
liftover=$yc/software/liftover/liftOver
chainfile_hg1819=$yc/software/liftover/ucsc/hg18ToHg19.over.chain.gz
chainfile_hg1938=$yc/software/liftover/ucsc/hg19ToHg38.over.chain.gz

#-------------------------------------------------------------------------#
#                         Trinn1: create intial set
#-------------------------------------------------------------------------#
run_flag='N'
if [ $run_flag = "Y" ]; then
cut -f7-8,12-13,22 $gwasc | awk -F '\t' '{print $3,$4,$5,$1,$2}' OFS='\t' | tr -dc '\0-\177' | sed 's/&beta;/Beta /g' > $gwasc.txt

#filter irregular snps with empty chr
cut -f1-3 $gwasc.txt | awk -F '\t' '$1==""' | grep -v ":" | grep -i -v rs | grep -i -v chr | cut -f3 > $gwasc.tmp
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $gwasc.tmp $gwasc.txt > $gwasc.tmp2
mv $gwasc.tmp2 $gwasc.txt

#clean snps
cut -f1-3 $gwasc.txt | sed 's/\./:/' | sed 's/\s+//g' | sed 's/ ; /;/g' | sed 's/; /;/g' | sed 's/ ;/;/g' | sed 's/ , /;/g' | sed 's/, /;/g' | sed 's/ ,/;/g' | sed 's/,/;/g' | sed 's/ \/ /,/g' | sed 's/ \//,/g' | sed 's/\/ /,/g' | sed 's/\//,/g' | sed 's/ x /x/gi' > $gwasc.tmp1
cut -f4-5 $gwasc.txt > $gwasc.tmp2
for((i=0; i<=22; i++)); do
    sed "s/chr${i}: /chr${i}:/gi" $gwasc.tmp1 > $gwasc.tmp3
    sed "s/chr${i}_/chr${i}:/gi" $gwasc.tmp3 > $gwasc.tmp1
done
sed 's/chr://gi' $gwasc.tmp1 | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/del-//gi' | sed 's/		.*psy_/		/gi' | sed 's/hg18_/hg18/gi' > $gwasc.tmp3
paste -d'	' $gwasc.tmp3 $gwasc.tmp2 > $gwasc.tmp

#process entries with empty chr/pos
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/' > $gwasc.tmp1
rm -f $gwasc.tmp10

cut -f3- $gwasc.tmp1 | while read line; do
    study=`echo "$line" | cut -f2-3`
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | while read line; do
        snp=`echo $line | cut -d'_' -f1 | cut -d'-' -f1 | cut -d':' -f1-2 | sed 's/[a-zA-Z]$//'`
        if [ `echo $snp | grep ':' | wc -l` -gt 0 ]; then
            chr=`echo $snp | cut -d':' -f1`
            pos=`echo $snp | cut -d':' -f2`
            echo $chr'	'$pos'	'$snp'	'"$study" >> $gwasc.tmp10
        else
            echo '-	-	'$snp'	'"$study" >> $gwasc.tmp10
        fi
    done 
done
mv $gwasc.tmp10 $gwasc.tmp1

#process entries with single snp with chr/pos
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/ && $3!~/x/' $gwasc.tmp | tail -n +2 | cut -f1-5 > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f4-5 $gwasc.tmp2 > $gwasc.tmp23
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 > $gwasc.tmp2

#process entries with multiple snps with chr/pos
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/ || $3~/x/)' $gwasc.tmp > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp31
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp32
    echo "$line" | cut -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | sed 's/[a-zA-Z]$//' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
    study=`echo "$line" | cut -f4-5`
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
    if [ $n3 -eq $n1 ] && [ `echo "$line" | cut -f3 | grep 'x' | wc -l` -gt 0 ]; then
        rm -f $gwasc.tmp34
        for ((i=1; i<=$n3; i++)); do
            echo "$study" >> $gwasc.tmp34
        done
    elif [ $n1 -eq 1 ] && [ $n3 -gt $n1 ] && [ `echo "$line" | cut -f3 | grep ',' | wc -l` -gt 0 ]; then
        for((i=2; i<=$n3; i++)); do
            head -n1 $gwasc.tmp31 >> $gwasc.tmp31
            head -n1 $gwasc.tmp32 >> $gwasc.tmp32
        done
        rm -f $gwasc.tmp34
        for ((i=1; i<=$n3; i++)); do
            echo "$study" >> $gwasc.tmp34
        done
    else
        paste -d':' $gwasc.tmp31 $gwasc.tmp32 > $gwasc.tmp34
        for((i=1; i<=$n3; i++)); do
            echo '-' >> $gwasc.tmp31
            echo '-' >> $gwasc.tmp32
        done
        cat $gwasc.tmp33 >> $gwasc.tmp34
        mv $gwasc.tmp34 $gwasc.tmp33
        rm -f $gwasc.tmp34
        for ((i=1; i<=$((n1+n3)); i++)); do
            echo "$study" >> $gwasc.tmp34
        done
    fi
    paste -d'	' $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp33 $gwasc.tmp34 >> $gwasc.tmp30
    rm -f $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp33 $gwasc.tmp34
done
mv $gwasc.tmp30 $gwasc.tmp3

echo 'CHR	POS	SNP	STUDY	TRAIT' > $outfile
cat $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3 | sed 's/^23	/X	/' | sed 's/	23:/	X:/' | sed 's/^24	/Y	/' | sed 's/	24:/	Y:/' | sort | uniq >> $outfile
rm -f $gwasc.tmp*

#update chr and pos of some empty entries according to available entries
tail -n +2 $outfile | grep -v ^- > $outfile.tmp
cut -f1-3 $outfile.tmp | sort -s -k3,3 | uniq > $outfile.tmp1
grep ^- $outfile | sort -s -k3,3 > $outfile.tmp2
join -1 3 -2 3 -a 1 -t '	' $outfile.tmp2 $outfile.tmp1 | awk -F '\t' '{if (NF>5) print $6,$7,$1,$4,$5; else print $2,$3,$1,$4,$5}' OFS='\t' >> $outfile.tmp

#compile whole set (OUTPUT: CHR POS SNP STUDY TRAIT)
head -n1 $outfile > $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1~/^[0-9]+$/' | sort -n -k1,1 -k2,2 >> $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1!~/^[0-9]+$/' | sort -k1,1 >> $outfile.tmp2
grep ^- $outfile.tmp >> $outfile.tmp2
mv $outfile.tmp2 $outfile
if [ `cut -f1-3 $outfile | grep : | sed 's/hg18//g' | awk '$1":"$2!=$3' | wc -l` -gt 0 ]; then
    echo "NOTE: SNPs with inconsistent chr/pos" 
    cut -f1-3 $outfile | grep : | sed 's/hg18//g' | awk '$1":"$2!=$3'
    exit 1
fi
rm -f $outfile.tmp*

#-------------------------------------------------------------------------#
#Trinn2: use build37_snp_chr_pos to find chr:pos for some snps ready for
#Batch 2
#-------------------------------------------------------------------------#
#snps with only rs
grep ^- $outfile | sort -s -k3,3 > $outfile.tmp

#try to convert rs to chr:pos (build37/hg19)
if [ ! -f $(dirname $outfile)/build37_SNPChrPos_map.txt ]; then
    zcat $build37_snpchrpos_ref | cut -f1-3 | awk 'NF==3' | sed 's/^/rs/' > $(dirname $outfile)/build37_SNPChrPos_map.txt
fi
awk -F '\t' 'NR==FNR {D[$3]++; next} ($1 in D)' $outfile.tmp $(dirname $outfile)/build37_SNPChrPos_map.txt | sort -s -k1,1 > $(dirname $outfile)/build37_SNPChrPos_map.tmp
join -1 3 -2 1 -a 1 -t '	' $outfile.tmp $(dirname $outfile)/build37_SNPChrPos_map.tmp | awk -F '\t' '{if(NF==7) print $6,"hg19"($7+1),$1,$4,$5; else print $2,$3,$1,$4,$5}' OFS='\t' | sort -s -k3,3 > $outfile.tmp2
#-------------------------------------------------------------------------#

#rm -f $outfile.tmp3
#grep ^- $outfile.tmp2 | while read line; do
#    rs=`echo "$line" | cut -f3`
#    echo "" >> $outfile.tmp3
#    echo $rs >> $outfile.tmp3
#    if [ ! -s index.html?term=$rs ]; then
#        wget https://www.ncbi.nlm.nih.gov/snp/?term=$rs
#    fi
#    snp=$rs
#    if [ -s index.html?term=$rs ]; then
#        if [ `grep Chromosome: index.html?term=$rs | wc -l` -eq 0 ]; then
#            continue
#        fi
#        if [ `grep Chromosome: index.html?term=$rs | grep 'no mapping' | wc -l` -gt 0 ]; then
#            continue
#        fi
#        if [ `grep 'Homo sapiens' index.html?term=$rs | wc -l` -gt 0 ]; then
#            if [ `grep 'Homo sapiens' index.html?term=$rs | grep "${rs}.*has merged into" | wc -l` -gt 0 ]; then
#                snp=`grep 'Homo sapiens' index.html?term=$rs | grep "${rs}.*has merged into" | sed 's/^.*has merged into//' | cut -d'>' -f2 | cut -d'<' -f1`
#            elif [ `grep 'Homo sapiens' index.html?term=$rs | grep -v "has merged into.*$rs" | wc -l` -gt 0 ]; then
#                snp=$rs
#            fi
#        fi
#        if [ `grep Chromosome: index.html?term=$rs | grep 'NT_' | wc -l` -gt 0 ]; then
#            chr=`grep Chromosome: index.html?term=$rs | cut -d'>' -f7 | cut -d: -f1 | head -n1`
#            pos38=`grep Chromosome: index.html?term=$rs | cut -d'>' -f7 | cut -d: -f2 | head -n1`
#        else
#            chr=`grep Chromosome: index.html?term=$rs | cut -d'>' -f7 | cut -d: -f1 | head -n1`
#            pos38=`grep Chromosome: index.html?term=$rs | cut -d'>' -f7 | cut -d: -f2 | head -n1`
#            pos37=`grep GRCh38 index.html?term=$rs | cut -d'>' -f2 | cut -d: -f2 | head -n1`
#        fi
#        if [ `echo $chr | grep ^MT | wc -l` -gt 0 ]; then
#            chr=M
#        fi
#        if [ `echo $chr | grep ^N | wc -l` -gt 0 ]; then
#            chr=6
#        fi
#        echo $snp $chr $pos38 >> $outfile.tmp3
#    fi
#done

#-------------------------------------------------------------------------#
#                    Trinn3: snps with only rs - Batch 1
#-------------------------------------------------------------------------#
#create RS liftover map
if [ ! -f $(dirname $outfile)/RS_liftover_map.txt ]; then
    zcat $yc/software/liftover/ncbi/SNPHistory.bcp.gz | grep -i -v activ | cut -f1 > $(dirname $outfile)/SNPHistory.txt
    zcat $yc/software/liftover/ncbi/RsMergeArch.bcp.gz | awk -F '\t' 'NF>=7' | cut -f1-2,7 | awk 'NF==3' | grep -v [a-zA-Z] | grep -v ':' | grep -v '-' | grep -v '\.' | grep -v ' ' > $(dirname $outfile)/RsMergeArch.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/SNPHistory.txt $(dirname $outfile)/RsMergeArch.txt | awk '$1!=$3 {print $1,$3}' OFS='\t' > $(dirname $outfile)/RS_liftover_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/SNPHistory.txt $(dirname $outfile)/RsMergeArch.txt | awk '$2!=$3 {print $2,$3}' OFS='\t' >> $(dirname $outfile)/RS_liftover_map.txt
    sort -n -k2,2 $(dirname $outfile)/RS_liftover_map.txt | uniq | sed 's/^/rs/' | sed 's/	/	rs/' > $(dirname $outfile)/RS_liftover_map.tmp
    mv $(dirname $outfile)/RS_liftover_map.tmp $(dirname $outfile)/RS_liftover_map.txt
fi

#snps with only rs info
grep ^- $outfile.tmp2 > $outfile.tmp3

#update rs with only one merging record (and leave the potential ones with multiple to be checked in the end without chr:pos)
awk -F '\t' 'NR==FNR {D[$3]++; next} ($1 in D)' $outfile.tmp3 $(dirname $outfile)/RS_liftover_map.txt | sort -s -k1,1 | awk '{print $2,$1}' | uniq -c -f1 | awk '$1==1 {print $3,$2}' OFS='\t' > $(dirname $outfile)/RS_liftover_map.tmp
join -1 3 -2 1 -a 1 -t '	' $outfile.tmp3 $(dirname $outfile)/RS_liftover_map.tmp | awk -F '\t' '{if(NF==6) print $2,$3,$6,$4,$5; else print $2,$3,$1,$4,$5}' OFS='\t' | sort -s -k3,3 > $outfile.tmp4
rm -f $(dirname $outfile)/RS_liftover_map.tmp

#make hg38 chr:pos snp map
mkdir -p $(dirname $outfile)/hg38
if [ ! -f $(dirname $outfile)/hg38/hg38_snps.csv ]; then
    zcat $hg38_chrpossnp_ref | cut -f2,3,5 | sed 's/^chr//' | awk -F '\t' '{print $1":"($2+1),$3}' > $(dirname $outfile)/hg38/hg38_snps.csv
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
        echo "chr"$i
        grep "^$chr:" $(dirname $outfile)/hg38/hg38_snps.csv | sed 's/ /	/' > $(dirname $outfile)/hg38/hg38_snps.chr$i.csv
    done
fi

#make hg19 chr:pos snp map
mkdir -p $(dirname $outfile)/hg19
if [ ! -f $(dirname $outfile)/hg19/hg19_snps.csv ]; then
    zcat $hg19_chrpossnp_ref | cut -f2,3,5 | sed 's/^chr//' | awk -F '\t' '{print $1":"($2+1),$3}' > $(dirname $outfile)/hg19/hg19_snps.csv
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
        echo "chr"$i
        grep "^$chr:" $(dirname $outfile)/hg19/hg19_snps.csv | sed 's/ /	/' > $(dirname $outfile)/hg19/hg19_snps.chr$i.csv
    done
fi

#search chr:pos for those snps with only rs info
mkdir -p $(dirname $outfile)/rs
rm -f $(dirname $outfile)/withdrawn_and_unfound_rs.txt
rm -f $(dirname $outfile)/hg38_snp_map.txt
rm -f $(dirname $outfile)/hg19_snp_map.txt
grep ^- $outfile.tmp4 | cut -f3 | sort -s | uniq | while read rs; do
    #retrieve dbSNP - Short Genetic Variations
    if [ ! -f $(dirname $outfile)/rs/$rs ]; then
        wget https://www.ncbi.nlm.nih.gov/snp/$rs -O $(dirname $outfile)/rs/$rs
    fi

    #remove withdrawn and unfound snps
    if [ `cat $(dirname $outfile)/rs/$rs | wc -l` -eq 0 ]; then
        echo $rs"	unfound" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        continue
    fi
    if [ `grep "${rs} was withdrawn" $(dirname $outfile)/rs/$rs | wc -l` -gt 0 ]; then
        echo $rs"	withdrawn" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        continue
    fi

    #extract chr info
    chr=''
    chr=`grep chr $(dirname $outfile)/rs/$rs | grep : | grep GRCh | head -n1 | cut -d'>' -f2 | cut -d: -f1 | sed 's/chr//' | sed 's/MT/M/'`
    if [ "$chr" = "" ]; then
        chr=`grep chr $(dirname $outfile)/rs/$rs | grep : | grep '()' | head -n1 | cut -d'>' -f2 | cut -d: -f1 | sed 's/chr//' | sed 's/MT/M/'`
    fi
    if [ "$chr" = "" ]; then
        if [ `grep 'N.*_' $(dirname $outfile)/rs/$rs | grep : | grep GRCh | wc -l` -gt 0 ]; then
            chr=`grep 'chr ' $(dirname $outfile)/rs/$rs | grep -v \>  | sed 's/^.*chr//' | head -n1 | awk '{print $1}'`
            pos38=`grep "$rs$" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f1 | cut -d: -f2`
            if [ "$pos38" != "" ]; then
                echo $rs'	'$chr'	'$pos38 >> $(dirname $outfile)/hg38_snp_map.txt
            fi
            pos19=`grep "$rs$" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f1 | cut -d: -f2`
            if [ "$pos19" != "" ]; then
                echo $rs'	'$chr'	'$pos19 >> $(dirname $outfile)/hg19_snp_map.txt
            fi
            if [ "$pos38" != "" ] || [ "$pos19" != "" ]; then
                continue
            fi
            pos=`grep 'N.*_' $(dirname $outfile)/rs/$rs | grep : | grep GRCh38 | cut -d: -f2 | cut -d' ' -f1`
            if [ "$pos" != "" ]; then
                echo $rs'	'$chr'	'$pos >> $(dirname $outfile)/hg38_snp_map.txt
            fi
            pos=`grep 'N.*_' $(dirname $outfile)/rs/$rs | grep : | grep GRCh37 | cut -d: -f2 | cut -d' ' -f1`
            if [ "$pos" != "" ]; then
                echo $rs'	'$chr'	'$pos >> $(dirname $outfile)/hg19_snp_map.txt
            fi
            continue
        fi
    fi

    hg38id=''
    hg19id=''
    if [ "$chr" != "" ]; then
        hg38id=`grep "$rs$" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f1`
    else
        for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
            hg38id=`grep "$rs$" $(dirname $outfile)/hg38/hg38_snps.chr$i.csv | cut -f1`
            if [ "$hg38id" != "" ]; then
                break
            fi
        done
    fi
    if [ "$hg38id" != "" ]; then
        chr=`echo $hg38id | cut -d: -f1`
        pos=`echo $hg38id | cut -d: -f2`
        echo $rs'	'$chr'	'$pos >> $(dirname $outfile)/hg38_snp_map.txt
    fi
    if [ "$chr" != "" ]; then
        hg19id=`grep "$rs$" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f1`
    else
        for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
            hg19id=`grep "$rs$" $(dirname $outfile)/hg19/hg19_snps.chr$i.csv | cut -f1`
            if [ "$hg19id" != "" ]; then
                break
            fi
        done
    fi
    if [ "$hg19id" != "" ]; then
        chr=`echo $hg19id | cut -d: -f1`
        pos=`echo $hg19id | cut -d: -f2`
        echo $rs'	'$chr'	'$pos >> $(dirname $outfile)/hg19_snp_map.txt
    fi
done

#remove withdrawn and unfound rs
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/withdrawn_and_unfound_rs.txt $outfile.tmp4 > $outfile.tmp5
fi
#update snps with hg19 chr:pos
join -1 3 -2 1 -a 1 -t '	' $outfile.tmp5 $(dirname $outfile)/hg19_snp_map.txt | awk -F '\t' '{if (NF==7) print $6":"$7,$6,$7,$1,$4,$5,$6":"$7,"-"; else print $2":"$3,$2,$3,$1,$4,$5,"-","-"}' OFS='\t' > $outfile.tmp6
#update snps with hg38 chr:pos (OUTPUT: CHR:POS CHR POS RS STUDY TRAIT Hg19 Hg38)
join -1 4 -2 1 -a 1 -t '	' $outfile.tmp6 $(dirname $outfile)/hg38_snp_map.txt | awk -F '\t' '{if (NF==10) print $9":"$10,$9,$10,$1,$5,$6,$7,$9":"$10; else print $2,$3,$4,$1,$5,$6,$7,$8}' OFS='\t' > $outfile.tmp7

if [ `awk -F '\t' '$7=="-" && $8=="-" {print $4}' $outfile.tmp7 | wc -l` -gt 0 ]; then
    echo "snps without retrieved hg19 and hg38 chr:pos"
    awk -F '\t' '$7=="-" && $8=="-" {print $4}' $outfile.tmp7
fi
n_snps=`grep ^- $outfile.tmp4 | cut -f3 | wc -l`
n_uniq_snps=`grep ^- $outfile.tmp4 | cut -f3 | sort | uniq | wc -l`
n_snps_remain=`cut -f4 $outfile.tmp7 | wc -l`
n_uniq_snps_remain=`cut -f4 $outfile.tmp7 | sort | uniq | wc -l`
n_withdraw_unfound=`cat $(dirname $outfile)/withdrawn_and_unfound_rs.txt | wc -l`
n_hg19_map=`cat hg19_snp_map.txt | wc -l`
n_hg38_map=`cat hg38_snp_map.txt | wc -l`
n_hg19_lack=`awk -F '\t' '$7=="-" {print $8}' $outfile.tmp7 | wc -l`
n_hg38_lack=`awk -F '\t' '$8=="-" {print $7}' $outfile.tmp7 | wc -l`
n_hg1938_lack=`awk -F '\t' '$7=="-" && $8=="-" {print $4}' $outfile.tmp7 | wc -l`
echo "ROUND 1"
echo "snps:$n_snps(uniq:$n_uniq_snps) withdraw/unfound:$((n_snps-n_snps_remain))(uniq:$n_withdraw_unfound) hg19_map(uniq:$n_hg19_map) hg38_map(uniq:$n_hg38_map)"
echo "snps_remain:$n_snps_remain(uniq:$n_uniq_snps_remain) hg19_lack:$n_hg19_lack(uniq:$((n_uniq_snps-n_withdraw_unfound-n_hg19_map))) hg38_lack:$n_hg38_lack(uniq:$((n_uniq_snps-n_withdraw_unfound-n_hg38_map))) hg19_hg38_lack:$n_hg1938_lack"

run_flag_2='N'
if [ $run_flag_2 = "Y" ]; then
#-------------------------------------------------------------------------#
#                    Trinn4: snps with chr:pos - Batch 2
#-------------------------------------------------------------------------#
#snps with chr:pos
tail -n +2 $outfile | grep -v ^- | awk '$3!~/hg18/' | awk -F '\t' '{print $1":"$2,$1,$2,$3,$4,$5,"-",$1":"$2}' OFS='\t' > $outfile.2.tmp
tail -n +2 $outfile | grep -v ^- | awk '$3~/hg18/' | sed 's/hg18//' | awk -F '\t' '{print $1":"$2,$1,$2,$3,$4,$5,"-","-"}' OFS='\t' >> $outfile.2.tmp
#snps with only rs and updated chr:pos from build37/hg19
grep -v ^- $outfile.tmp2 | sed 's/hg19//g' | awk -F '\t' '{print "hg19"$1":"$2,$1,$2,$3,$4,$5,$1":"$2,"-"}' OFS='\t' | sed 's/:hg19/:/' >> $outfile.2.tmp
sort $outfile.2.tmp | uniq | sort -s -k1,1 > $outfile.2.tmp2
mv $outfile.2.tmp2 $outfile.2.tmp

#run UCSC liftOver for positional map (hg18->hg19)
tail -n +2 $outfile | cut -f1-3 | grep -v ^- | grep hg18 | sed 's/hg18//' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}.hg18.txt
if [ -s ${outfile%.*}.hg18.txt ]; then
    $liftover ${outfile%.*}.hg18.txt $chainfile_hg1819 ${outfile%.*}_liftover_hg18_out_hg19.bed ${outfile%.*}_liftover_unlifted_hg18.bed
    cat ${outfile%.*}_liftover_hg18_out_hg19.bed | awk '{print $4,$2}' OFS='\t' >  ${outfile%.*}_liftover_hg19_map.txt
    cat ${outfile%.*}_liftover_unlifted_hg18.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg18_snps.txt

    #delete unlifted hg18 snps 
    if [ -s ${outfile%.*}_liftover_hg18_unlifted_snps.txt ]; then
        awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' ${outfile%.*}_liftover_hg18_unlifted_snps.txt $outfile.2.tmp > $outfile.2.tmp2
        mv $outfile.2.tmp2 $outfile.2.tmp
    fi

    #update snps with lifted hg19 positions
    if [ -s ${outfile%.*}_liftover_hg19_map.txt ]; then
        sort ${outfile%.*}_liftover_hg19_map.txt | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg19_map.tmp
        join -1 1 -2 1 -a 1 -t '	' $outfile.2.tmp ${outfile%.*}_liftover_hg19_map.tmp | awk -F '\t' '{if(NF==9) print "hg19"$2":"$9,$2,$9,$4,$5,$6,$2":"$9,$8; else print $0}' OFS='\t' > $outfile.2.tmp2
        mv $outfile.2.tmp2 $outfile.2.tmp
        rm -f ${outfile%.*}_liftover_hg19_map.tmp
    fi
fi

#snps with chr:pos
awk '$1!~/hg19/' $outfile.2.tmp > $outfile.2.tmp2
awk '$1~/hg19/' $outfile.2.tmp | sed 's/hg19//' >> $outfile.2.tmp2
sort $outfile.2.tmp2 | uniq | sort -s -k1,1 > $outfile.2.tmp3
mv $outfile.2.tmp3 $outfile.2.tmp2

#run UCSC liftOver for positional map (hg19->hg38)
cut -f1-3 $outfile.2.tmp | grep hg19 | sed 's/hg19//' | awk '{print "chr"$2,$3,($3+1),$1}' OFS='\t' > ${outfile%.*}.hg19.txt
if [ -s ${outfile%.*}.hg19.txt ]; then
    $liftover ${outfile%.*}.hg19.txt $chainfile_hg1938 ${outfile%.*}_liftover_hg19_out_hg38.bed ${outfile%.*}_liftover_unlifted_hg19.bed
    cat ${outfile%.*}_liftover_hg19_out_hg38.bed | awk '{print $4,$2}' OFS='\t' >  ${outfile%.*}_liftover_hg38_map.txt
    cat ${outfile%.*}_liftover_unlifted_hg19.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg19_snps.txt

    #delete unlifted hg19 snps 
    if [ -s ${outfile%.*}_liftover_unlifted_hg19_snps.txt ]; then
        rm -f ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt
        cat ${outfile%.*}_liftover_unlifted_hg19_snps.txt | while read snp; do
            chr=${snp%:*}
            rs=`grep "^$snp	" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f2`
            if [ "$rs" = "" ]; then
                echo $snp >> ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt
                continue 
            fi
            #retrieve dbSNP - Short Genetic Variations
            #if [ ! -f $(dirname $outfile)/rs/$rs ]; then
            #    wget https://www.ncbi.nlm.nih.gov/snp/$rs -O $(dirname $outfile)/rs/$rs
            #fi
            #remove withdrawn and unfound snps
            #if [ `cat $(dirname $outfile)/rs/$rs | wc -l` -eq 0 ]; then
            #    echo $snp >> ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt
            #    continue
            #fi
            #if [ `grep "${rs} was withdrawn" $(dirname $outfile)/rs/$rs | wc -l` -gt 0 ]; then
            #    echo $snp >> ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt       
            #    continue
            #elif [ `grep chr $(dirname $outfile)/rs/$rs | grep : | grep GRCh | wc -l` -gt 0 ]; then
            #    continue
            #elif [ `grep chr $(dirname $outfile)/rs/$rs | grep : | grep '()' | wc -l` -gt 0 ]; then
            #    continue
            #elif [ `grep 'N.*_' $(dirname $outfile)/rs/$rs | grep : | grep GRCh | wc -l` -gt 0 ]; then
            #    continue
            #else
            #    echo $snp >> ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt
            #fi
        done
        if [ -s  ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt ]; then
            awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' ${outfile%.*}_liftover_unlifted_hg19_snps_2.txt $outfile.2.tmp2 > $outfile.2.tmp3
            mv $outfile.2.tmp3 $outfile.2.tmp2
        fi
    fi
    
    n_snps=`cat $outfile.2.tmp2 | wc -l`
    n_hg19=`awk -F '\t' '$7==$1' $outfile.2.tmp2 | wc -l`
    n_hg38_lack=`awk -F '\t' '$8=="-"' $outfile.2.tmp2 | wc -l`
    n_hg38=`awk -F '\t' '$8==$1' $outfile.2.tmp2 | wc -l`
    n_hg19_lack=`awk -F '\t' '$7=="-"' $outfile.2.tmp2 | wc -l`
    echo "snps:$n_snps hg19:$n_hg19($n_hg38_lack) hg38:$n_hg38($n_hg19_lack)"

    #update snps with lifted hg38 positions (OUTPUT: CHR:POS CHR POS SNP STUDY TRAIT Hg19 Hg38)
    if [ -s ${outfile%.*}_liftover_hg38_map.txt ]; then
        sort ${outfile%.*}_liftover_hg38_map.txt | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg38_map.tmp
        join -1 1 -2 1 -a 1 -t '	' $outfile.2.tmp2 ${outfile%.*}_liftover_hg38_map.tmp | awk -F '\t' '{if(NF==9) print $2":"$9,$2,$9,$4,$5,$6,$7,$2":"$9; else print $0}' OFS='\t' > $outfile.2.tmp3
        rm -f ${outfile%.*}_liftover_hg38_map.tmp
    fi
    n_hg1938=`awk -F '\t' '$8==$7 {print $8}' $outfile.2.tmp3 | wc -l`
    n_hg3819=`sed 's/:/ /' ${outfile%.*}_liftover_hg38_map.txt | awk '$2==$3' | sort | uniq | wc -l`
    n_hg19_lack=`awk -F '\t' '$7=="-" {print $8}' $outfile.2.tmp3 | wc -l`
    n_hg38_lack=`awk -F '\t' '$8=="-" {print $7}' $outfile.2.tmp3 | wc -l`
    echo "ROUND 2"
    n_snps_remain=`cut -f4 $outfile.2.tmp3 | wc -l`
    echo "hg19_hg39_shared_pos:$n_hg1938(uniq:$n_hg3819)"
    echo "hg19_lack:$n_hg19_lack hg38_lack:$n_hg38_lack"
    awk -F '\t' '$8=="-" {print "hg19:"$7}' $outfile.2.tmp3
fi
fi
#rm -f $outfile*.tmp*
#-------------------------------------------------------------------------#
