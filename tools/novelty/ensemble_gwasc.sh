#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to ensemble NHGRI-EBI GWAS Catalog to make a clean file
# for novelty checking of loci WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage: sh ensemble_gwasc.sh gwasc outfile"
  echo "Arguments: gwasc - gwascatalog"
  echo "           outfile - output file"
  echo "Example: sh ensemble_gwasc.sh gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv gwas_catalog_v1.0-associations_e100_r2021-02-25.csv"
  exit 0
fi
#-------------------------------------------------------------------------#

gwasc=$1
outfile=$2

build37_snpchrpos_ref=$yc/software/liftover/ncbi/b132_SNPChrPosOnRef_37_1.bcp.gz
build38_snpchrpos_ref=$yc/software/liftover/ncbi/b151_SNPChrPosOnRef_108.bcp.gz
hg19_chrpossnp_ref=$yc/software/liftover/ucsc/hg19/snp151.txt.gz
hg38_chrpossnp_ref=$yc/software/liftover/ucsc/hg38/snp151.txt.gz
liftover=$yc/software/liftover/liftOver
chainfile_hg1819=$yc/software/liftover/ucsc/hg18ToHg19.over.chain.gz
chainfile_hg1938=$yc/software/liftover/ucsc/hg19ToHg38.over.chain.gz
chainfile_hg3819=$yc/software/liftover/ucsc/hg38ToHg19.over.chain.gz

#-------------------------------------------------------------------------#
#                    Trinn1: create maps
#-------------------------------------------------------------------------#
mkdir -p $(dirname $outfile)/map
#create build37 chr:pos snp map
if [ ! -f $(dirname $outfile)/map/build37_SNPChrPos_map.txt ]; then
    zcat $build37_snpchrpos_ref | cut -f1-3 | awk 'NF==3' | sed 's/^/rs/' | awk '{print $1,$2,($3+1)}' OFS='\t' > $(dirname $outfile)/map/build37_SNPChrPos_map.txt
fi
#create build38 chr:pos snp map
mkdir -p $(dirname $outfile)/b38
if [ ! -f $(dirname $outfile)/b38/b38_snps.csv ]; then
    zcat $build38_snpchrpos_ref | cut -f1-3 | awk 'NF==3' | sed 's/^/rs/' | sed 's/MT/M/' | awk '{print $2,($3+1),$1}' OFS='\t' > $(dirname $outfile)/b38/b38_snps.csv
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
        echo "chr"$i
        grep "^$i	" $(dirname $outfile)/b38/b38_snps.csv | awk -F '\t' '{print $1":"$2,$1,$2,$3}' OFS='\t' > $(dirname $outfile)/b38/build38_snps.chr$i.csv
    done
fi

#create RS liftover map
if [ ! -f $(dirname $outfile)/map/RS_liftover_map.txt ]; then
    zcat $yc/software/liftover/ncbi/SNPHistory.bcp.gz | grep -i -v activ | cut -f1 > $(dirname $outfile)/map/SNPHistory.txt
    zcat $yc/software/liftover/ncbi/RsMergeArch.bcp.gz | awk -F '\t' 'NF>=7' | cut -f1-2,7 | awk 'NF==3' | grep -v [a-zA-Z] | grep -v ':' | grep -v '-' | grep -v '\.' | grep -v ' ' > $(dirname $outfile)/map/RsMergeArch.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/map/SNPHistory.txt $(dirname $outfile)/map/RsMergeArch.txt | awk '$1!=$3 {print $1,$3}' OFS='\t' > $(dirname $outfile)/map/RS_liftover_map.txt
    awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/map/SNPHistory.txt $(dirname $outfile)/map/RsMergeArch.txt | awk '$2!=$3 {print $2,$3}' OFS='\t' >> $(dirname $outfile)/map/RS_liftover_map.txt
    sort -n -k2,2 $(dirname $outfile)/map/RS_liftover_map.txt | uniq | sed 's/^/rs/' | sed 's/	/	rs/' > $(dirname $outfile)/map/RS_liftover_map.tmp
    mv $(dirname $outfile)/map/RS_liftover_map.tmp $(dirname $outfile)/map/RS_liftover_map.txt
fi

#make hg38 chr:pos snp map
mkdir -p $(dirname $outfile)/hg38
if [ ! -f $(dirname $outfile)/hg38/hg38_snps.csv ]; then
    zcat $hg38_chrpossnp_ref | cut -f2-5,11-12,17 | sed 's/^chr//' > $(dirname $outfile)/hg38/hg38_snps.csv
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
        echo "chr"$i
        #grep "^$i	" $(dirname $outfile)/hg38/hg38_snps.csv | awk -F '\t' '{if($2==$3) print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$1":"$2; else print $1":"($2+1),$1,$2,$3,$4,$5,$6,$7,$1":"$2}' OFS='\t' > $(dirname $outfile)/hg38/hg38_snps.chr$i.csv
        grep "^$i	" $(dirname $outfile)/hg38/hg38_snps.csv | awk -F '\t' '{print $1":"($2+1),$1,$2,$3,$4,$5,$6,$7,$1":"$2}' OFS='\t' > $(dirname $outfile)/hg38/hg38_snps.chr$i.csv
    done
fi

#make hg19 chr:pos snp map
mkdir -p $(dirname $outfile)/hg19
if [ ! -f $(dirname $outfile)/hg19/hg19_snps.csv ]; then
    zcat $hg19_chrpossnp_ref | cut -f2-5,11-12,17 | sed 's/^chr//' > $(dirname $outfile)/hg19/hg19_snps.csv
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
        echo "chr"$i
        #grep "^$i	" $(dirname $outfile)/hg19/hg19_snps.csv | awk -F '\t' '{if($2==$3) print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$1":"$2; else print $1":"($2+1),$1,$2,$3,$4,$5,$6,$7,$1":"$2}' OFS='\t' > $(dirname $outfile)/hg19/hg19_snps.chr$i.csv
        grep "^$i	" $(dirname $outfile)/hg19/hg19_snps.csv | awk -F '\t' '{print $1":"($2+1),$1,$2,$3,$4,$5,$6,$7,$1":"$2}' OFS='\t' > $(dirname $outfile)/hg19/hg19_snps.chr$i.csv
    done
fi
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#                         Trinn2: create initial set
#-------------------------------------------------------------------------#
run_flag='N'
if [ $run_flag = "Y" ]; then
echo "ROUND I compile initial set"
cut -f7-8,12-13,22 $gwasc | awk -F '\t' '{print $3,$4,$5,$1,$2}' OFS='\t' | tr -dc '\0-\177' | sed 's/&beta;/Beta /g' > $gwasc.txt

#filter irregular snps with empty chr
echo "filter irregular snps with empty chr"
cut -f1-3 $gwasc.txt | awk -F '\t' '$1==""' | grep -v ":" | grep -i -v rs | grep -i -v chr | cut -f3 > $gwasc.tmp
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $gwasc.tmp $gwasc.txt > $gwasc.tmp2
mv $gwasc.tmp2 $gwasc.txt

#mark: need to verify if orders of rs delimited with 'x' in line with chr:pos

#clean snps
echo "clean snps"
cut -f1-3 $gwasc.txt | sed 's/\./:/' | sed 's/\s+//g' | sed 's/ ; /;/g' | sed 's/; /;/g' | sed 's/ ;/;/g' | sed 's/ , /;/g' | sed 's/, /;/g' | sed 's/ ,/;/g' | sed 's/,/;/g' | sed 's/ \/ /,/g' | sed 's/ \//,/g' | sed 's/\/ /,/g' | sed 's/\//,/g' | sed 's/ x /x/gi' > $gwasc.tmp1
cut -f4-5 $gwasc.txt > $gwasc.tmp2
for((i=0; i<=22; i++)); do
    sed "s/chr${i}: /chr${i}:/gi" $gwasc.tmp1 > $gwasc.tmp3
    sed "s/chr${i}_/chr${i}:/gi" $gwasc.tmp3 > $gwasc.tmp1
done
sed 's/chr://gi' $gwasc.tmp1 | sed 's/chr//gi' | sed 's/che//gi' | sed 's/ch//gi' | sed 's/		.*psy_/		/gi' | sed 's/hg18_/hg18/gi' > $gwasc.tmp3
paste -d'	' $gwasc.tmp3 $gwasc.tmp2 > $gwasc.tmp

#process entries with empty chr:pos
echo "process entries with empty chr:pos"
awk -F '\t' '$1==""' $gwasc.tmp | awk -F '\t' '$3!~/*/ && $3!~/im/' > $gwasc.tmp1
rm -f $gwasc.tmp10

cut -f3- $gwasc.tmp1 | while read line; do
    study=`echo "$line" | cut -f2-3`
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | while read id; do
        snp=`echo $id | sed 's/del-//gi' | cut -d'_' -f1 | cut -d'-' -f1 | cut -d':' -f1-2 | sed 's/[a-zA-Z]$//'`
        if [ `echo $snp | grep ':' | wc -l` -gt 0 ]; then
            chr=`echo $snp | cut -d':' -f1`
            pos=`echo $snp | cut -d':' -f2`
            echo $chr'	'$pos'	'$snp'	'$id'	'"$study" >> $gwasc.tmp10
        else
            echo '-	-	'$snp'	'$id'	'"$study" >> $gwasc.tmp10
        fi
    done 
done
mv $gwasc.tmp10 $gwasc.tmp1

#process entries with single snp with chr:pos
echo "process entries with single snp with chr:pos"
awk -F '\t' '$1!="" && $3!~/;/ && $3!~/,/ && $3!~/x/' $gwasc.tmp | tail -n +2 | cut -f1-5 > $gwasc.tmp2
cut -f1-2 $gwasc.tmp2 > $gwasc.tmp21
cut -f3 $gwasc.tmp2 | sed 's/del-//gi'| sed 's/[a-zA-Z]$//' > $gwasc.tmp22
cut -f3 $gwasc.tmp2 > $gwasc.tmp23
cut -f4-5 $gwasc.tmp2 > $gwasc.tmp24
paste -d'	' $gwasc.tmp21 $gwasc.tmp22 $gwasc.tmp23 $gwasc.tmp24 > $gwasc.tmp2

#process entries with multiple snps with chr:pos
echo "process entries with multiple snps with chr:pos"
awk -F '\t' '$1!="" && ($3~/;/ || $3~/,/ || $3~/x/)' $gwasc.tmp > $gwasc.tmp3
rm -f $gwasc.tmp30
cat $gwasc.tmp3 | while read line; do
    echo "$line" | cut -f1 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp31
    echo "$line" | cut -f2 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' > $gwasc.tmp32
    echo "$line" | cut -f3 | sed 's/;/\n/g' | sed 's/,/\n/g' | sed 's/x/\n/g' | awk '$1~/:/ || $1~/rs/' > $gwasc.tmp33
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
        rm -f $gwasc.tmp31
        rm -f $gwasc.tmp32
        rm -f $gwasc.tmp34
        for((i=1; i<=$n3; i++)); do
            echo '-' >> $gwasc.tmp31
            echo '-' >> $gwasc.tmp32
            echo "$study" >> $gwasc.tmp34
        done
    fi
    sed 's/del-//gi' $gwasc.tmp33 | sed 's/[a-zA-Z]$//' > $gwasc.tmp35
    paste -d'	' $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp35 $gwasc.tmp33 $gwasc.tmp34 >> $gwasc.tmp30
    rm -f $gwasc.tmp31 $gwasc.tmp32 $gwasc.tmp35 $gwasc.tmp33 $gwasc.tmp34
done
mv $gwasc.tmp30 $gwasc.tmp3

cat $gwasc.tmp1 $gwasc.tmp2 $gwasc.tmp3 | sed 's/^23	/X	/' | sed 's/	23:/	X:/' | sed 's/^24	/Y	/' | sed 's/	24:/	Y:/' | sort | uniq > $outfile

#update some entries with only rs for chr:pos according to available entries
echo "update empty entries for chr:pos"
tail -n +2 $outfile | grep -v ^- > $outfile.tmp
cut -f1-3 $outfile.tmp | sort -s -k3,3 | uniq > $outfile.tmp1
grep ^- $outfile | sort -s -k3,3 > $outfile.tmp2
join -1 3 -2 3 -a 1 -t '	' $outfile.tmp2 $outfile.tmp1 | awk -F '\t' '{if (NF>6) print $7,$8,$1,$4,$5,$6; else print $2,$3,$1,$4,$5,$6}' OFS='\t' >> $outfile.tmp

#compile initial set (OUTPUT: CHR POS SNP ID STUDY TRAIT)
echo "compile initial set"
head -n1 $outfile > $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1~/^[0-9]+$/' | sort -n -k1,1 -k2,2 >> $outfile.tmp2
grep -v ^- $outfile.tmp | awk -F '\t' '$1!~/^[0-9]+$/' | sort -k1,1 >> $outfile.tmp2
grep ^- $outfile.tmp >> $outfile.tmp2
echo 'CHR	POS	SNP	ID	STUDY	TRAIT	FLAG' > $outfile
#Flag: 0-snp with only rs; 1-snp with chr:pos; 2-snp with both rs and chr:pos
awk -F '\t' '{if($1=="-") print $1,$2,$3,$4,$5,$6,"0"; else if($1!="-"&&$3~/:/) print $1,$2,$3,$4,$5,$6,"1"; else if($1!="-"&&$3~/rs/) print $1,$2,$3,$4,$5,$6,"2"; else print $1,$2,$3,$4,$5,$6,"x"}' OFS='\t' $outfile.tmp2 >> $outfile

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
fi

echo "ROUND II batch-wise process"
#-------------------------------------------------------------------------#
#                    Trinn3: snps with rs - Batch 1
#-------------------------------------------------------------------------#
run_flag='Y'
if [ $run_flag = "Y" ]; then
#search dbSNP chr:pos for those snps with rs info
mkdir -p $(dirname $outfile)/dbsnp
rm -f $(dirname $outfile)/withdrawn_and_unfound_rs.txt
rm -f $(dirname $outfile)/hg38_snp_map.txt
rm -f $(dirname $outfile)/hg19_snp_map.txt
rm -f $(dirname $outfile)/RS_merge_map.txt

awk -F '\t' '$7=="0" || $7=="2"' $outfile | cut -f3 | sort | uniq | while read rs; do
    #retrieve dbSNP - Short Genetic Variations
    if [ ! -f $(dirname $outfile)/dbsnp/$rs ]; then
        wget https://www.ncbi.nlm.nih.gov/snp/$rs -O $(dirname $outfile)/dbsnp/$rs
    fi

    #unfound snp
    if [ `cat $(dirname $outfile)/dbsnp/$rs | wc -l` -eq 0 ]; then
        echo $rs"	unfound" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        continue
    fi

    #withdrawn snp
    if [ `grep "${rs} was withdrawn" $(dirname $outfile)/dbsnp/$rs | wc -l` -gt 0 ]; then
        echo $rs"	withdrawn" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        continue
    fi

    #unsupported snp
    if [ `grep "no longer has any supporting" $(dirname $outfile)/dbsnp/$rs | wc -l` -gt 0 ] && [ `grep "former support of this RefSNP have moved to" $(dirname $outfile)/dbsnp/$rs | wc -l` -eq 0 ]; then
        echo $rs"	unsupported" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        continue
    fi

    snp=$rs
    #former support of observation
    if [ `grep "former support of this RefSNP have moved to" $(dirname $outfile)/dbsnp/$rs | wc -l` -gt 0 ]; then
        snp=`grep -A2 'former support of this RefSNP have moved to' $(dirname $outfile)/dbsnp/$rs | tail -n1 | cut -d'>' -f2 | cut -d'<' -f1`
        echo "$rs|$snp	formersupport" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        if [ ! -f $(dirname $outfile)/dbsnp/$snp ]; then
            wget https://www.ncbi.nlm.nih.gov/snp/$snp -O $(dirname $outfile)/dbsnp/$snp
        fi
        if [ -s $(dirname $outfile)/dbsnp/$snp ]; then  
            echo "$rs|$snp	formersupport" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        else
            echo $rs"	failed" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
        fi
    fi

    #merge into another
    if [ `grep "${rs} was merged into" $(dirname $outfile)/dbsnp/$rs | wc -l` -gt 0 ]; then
        snp=`grep -A3 "${rs} was merged into" dbsnp/$rs | tail -n1 | cut -d'>' -f2 | cut -d'<' -f1`
        if [ "$snp" != "" ]; then
            if [ ! -f $(dirname $outfile)/dbsnp/$snp ]; then
                wget https://www.ncbi.nlm.nih.gov/snp/$snp -O $(dirname $outfile)/dbsnp/$snp
            fi
            if [ -s $(dirname $outfile)/dbsnp/$snp ]; then  
                echo $rs'	'$snp >> $(dirname $outfile)/RS_merge_map.txt
            else
                echo $rs"	failed" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
            fi
        else
            snp=$rs
            echo $rs"	failed" >> $(dirname $outfile)/withdrawn_and_unfound_rs.txt
            continue
        fi
    fi

    #extract chr:pos info
    chr=`grep chr $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh | head -n1 | cut -d'>' -f2 | cut -d: -f1 | sed 's/chr//' | sed 's/MT/M/'`
    pos38=`grep chr $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh38 | head -n1 | cut -d'>' -f2 | cut -d: -f2 | cut -d' ' -f1 | cut -d'(' -f1 | cut -d'-' -f1`
    pos19=`grep chr $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh37 | head -n1 | cut -d'>' -f2 | cut -d: -f2 | cut -d' ' -f1 | cut -d'(' -f1 | cut -d'-' -f1`
    if [ "$chr" = "" ]; then
        chr=`grep chr $(dirname $outfile)/dbsnp/$snp | grep : | grep '()' | head -n1 | cut -d'>' -f2 | cut -d: -f1 | sed 's/chr//' | sed 's/MT/M/'`
    fi
    if [ "$chr" = "" ]; then
        if [ `grep 'N.*_' $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh | wc -l` -gt 0 ]; then
            chr=`grep 'chr ' $(dirname $outfile)/dbsnp/$snp | grep -v \>  | sed 's/^.*chr//' | head -n1 | awk '{print $1}'`
            pos38=`grep 'N.*_' $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh38 | cut -d: -f2 | cut -d' ' -f1`
            pos19=`grep 'N.*_' $(dirname $outfile)/dbsnp/$snp | grep : | grep GRCh37 | cut -d: -f2 | cut -d' ' -f1`
        fi
    fi
    if [ "$chr" != "" ]; then
        run_term_flag='N'
        if [ "$run_term_flag" = 'Y' ]; then
            if [ ! -s $(dirname $outfile)/dbsnp/term=$snp ]; then
                wget https://www.ncbi.nlm.nih.gov/snp/?term=$snp -O $(dirname $outfile)/dbsnp/term=$snp
            fi
            if [ -s $(dirname $outfile)/dbsnp/term=$snp ] && [ `grep Chromosome: $(dirname $outfile)/dbsnp/term=$snp | grep -v 'no mapping' | wc -l` -gt 0 ]; then
                if [ `grep 'Homo sapiens' $(dirname $outfile)/dbsnp/term=$snp | wc -l` -gt 0 ]; then
                    if [ `grep 'Homo sapiens' $(dirname $outfile)/dbsnp/term=$snp | grep "${snp}.*has merged into" | wc -l` -gt 0 ]; then
                        snp2=`grep 'Homo sapiens' $(dirname $outfile)/dbsnp/term=$snp | grep "${snp}.*has merged into" | sed 's/^.*has merged into//' | cut -d'>' -f2 | cut -d'<' -f1`
                    elif [ `grep 'Homo sapiens' $(dirname $outfile)/dbsnp/term=$snp | grep -v "has merged into.*$snp" | wc -l` -gt 0 ]; then
                        snp2=$snp
                    fi
                    if [ "$snp2" != "$snp" ] && [ "$snp" != "$rs" ]; then
                        echo "$rs|$snp merged into $snp2"
                        echo $snp'	'$snp2 >> $(dirname $outfile)/RS_liftover_map_2.txt
                        echo $rs'	'$snp2 >> $(dirname $outfile)/RS_liftover_map_2.txt
                    fi
                fi
            fi
            if [ `grep Chromosome: $(dirname $outfile)/dbsnp/term=$snp | grep 'no mapping' | wc -l` -eq 0 ]; then
                p38=`grep Chromosome: $(dirname $outfile)/dbsnp/term=$snp | cut -d'>' -f7 | cut -d: -f2 | head -n1`
                p19=`grep GRCh38 $(dirname $outfile)/dbsnp/term=$snp | cut -d'>' -f2 | cut -d: -f2 | head -n1`
                echo "$rs|$snp2"
                echo "pos38:$pos38	pos19:$pos19"
                echo "p38  :$p38	p19  :$p19"
                #if [ "$pos19" == "" ] && [ "$snp2" = "$rs" ] && [ "$chr" != "M" ]; then 
                #    if [ "$p19" != "" ] && [ "$p19" != "-1" ]; then
                #        pos19=$p19
                #    fi
                #fi
            fi
        fi
        if [ "$pos38" != "" ]; then 
            echo $rs'	'$chr'	'$pos38 >> $(dirname $outfile)/hg38_snp_map.txt
        fi
        if [ "$pos19" != "" ]; then 
            echo $rs'	'$chr'	'$pos19 >> $(dirname $outfile)/hg19_snp_map.txt
        fi
    fi
done
fi

run_flag='N'
if [ $run_flag = "Y" ]; then
#remove withdrawn and unfound rs
awk -F '\t' 'NR==FNR {D[$1]++; next} !($3 in D)' $(dirname $outfile)/withdrawn_and_unfound_rs.txt $outfile.tmp4 > $outfile.tmp5
#update snps with hg38 chr:pos
join -1 3 -2 1 -a 1 -t '	' $outfile.tmp5 $(dirname $outfile)/hg38_snp_map.txt | awk -F '\t' '{if (NF==7) print $6,$7,$1,$4,$5,"-",$6":"$7; else print $2,$3,$1,$4,$5,"-","-"}' OFS='\t' > $outfile.tmp6
#update snps with hg19 chr:pos (OUTPUT: CHR:POS CHR POS RS STUDY TRAIT Hg19 Hg38)
join -1 3 -2 1 -a 1 -t '	' $outfile.tmp6 $(dirname $outfile)/hg19_snp_map.txt | awk -F '\t' '{if (NF==9) print $8,$9,$1,$4,$5,$8":"$9,$7; else print $2,$3,$1,$4,$5,$6,$7}' OFS='\t' | awk -F '\t' '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7}' OFS='\t' | sort | uniq | sort -n -k2,2 -k3,3 > $outfile.tmp7

if [ `awk -F '\t' '$7=="-" && $8=="-" {print $4}' $outfile.tmp7 | wc -l` -gt 0 ]; then
    echo "snps_lack_both_hg19_and_hg38_chr:pos"
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
echo "n_snps:$n_snps(uniq:$n_uniq_snps) n_withdraw|unfound:$((n_snps-n_snps_remain))(uniq:$n_withdraw_unfound) n_hg19_map(uniq:$n_hg19_map) n_hg38_map(uniq:$n_hg38_map)"
echo "n_snps_remain:$n_snps_remain(uniq:$n_uniq_snps_remain) n_hg19_lack:$n_hg19_lack(uniq:$((n_uniq_snps-n_withdraw_unfound-n_hg19_map))) n_hg38_lack:$n_hg38_lack(uniq:$((n_uniq_snps-n_withdraw_unfound-n_hg38_map))) n_hg19_hg38_lack:$n_hg1938_lack"
echo "#-----------------------------------------------------------------------#"

#run UCSC liftOver for positional map (hg19->hg38)
awk -F '\t' '$7!="-" {print $7,$1}' OFS='\t' $outfile.tmp7 | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$3}' > ${outfile%.*}.hg19_0.txt
awk -F '\t' '$8!="-" {print $8,$1}' OFS='\t' $outfile.tmp7 | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$3}' > ${outfile%.*}.hg38_0.txt
sort -s -k1,1 $outfile.tmp7 >  $outfile.tmp8
if [ -s ${outfile%.*}.hg19_0.txt ]; then
    $liftover ${outfile%.*}.hg19_0.txt $chainfile_hg1938 ${outfile%.*}_liftover_hg19_out_hg38_0.bed ${outfile%.*}_liftover_unlifted_hg19_0.bed
    cat ${outfile%.*}_liftover_hg19_out_hg38_0.bed | awk '{print $4,$2}' OFS='\t' | sort | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg19_to_hg38_map_0.txt
    join -1 1 -2 1 -a 1 $outfile.tmp8 ${outfile%.*}_liftover_hg19_to_hg38_map_0.txt -t '	' | awk -F '\t' '{if(NF==9) print $1,$2,$3,$4,$5,$6,$7,$8,"-",$2":"$9; else print $1,$2,$3,$4,$5,$6,$7,$8,"-","-"}' OFS='\t' > $outfile.tmp9
    mv $outfile.tmp9 $outfile.tmp8
fi
if [ -s ${outfile%.*}.hg38_0.txt ]; then
    $liftover ${outfile%.*}.hg38_0.txt $chainfile_hg3819 ${outfile%.*}_liftover_hg38_out_hg19_0.bed ${outfile%.*}_liftover_unlifted_hg38_0.bed
    cat ${outfile%.*}_liftover_hg38_out_hg19_0.bed | awk '{print $4,$2}' OFS='\t' | sort | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg38_to_hg19_map_0.txt
    join -1 1 -2 1 -a 1 $outfile.tmp8 ${outfile%.*}_liftover_hg38_to_hg19_map_0.txt -t '	' | awk -F '\t' '{if(NF==11) print $1,$2,$3,$4,$5,$6,$7,$8,$2":"$11,$10; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' > $outfile.tmp9
    mv $outfile.tmp9 $outfile.tmp8
fi

awk -F '\t' '{if ($7==$9&&$8==$10) print $1,$2,$3,$4,$5,$6,$7,$8; else if($7=="-"&&$8!="-"&&$9!="-") print $1,$2,$3,$4,$5,$6,$9,$8; else if($7!="-") print $1,$2,$3,$4,$5,$6,$7,$10}' OFS='\t' $outfile.tmp8 | awk -F '\t' '{if($1==$7) print $1,$2,$3,$4,$5,$6,$7,$8,"hg19"; else print $1,$2,$3,$4,$5,$6,$7,$8,"hg38"}' OFS='\t' > $outfile.tmp9
fi

#-------------------------------------------------------------------------#
#                    Trinn4: snps with chr:pos - Batch 2
#-------------------------------------------------------------------------#
run_flag='N'
if [ $run_flag = "Y" ]; then
echo "Batch 2"
#snps with chr:pos
tail -n +2 $outfile | grep -v ^- | awk '$3!~/hg18/ && $3!~/hg19/' | awk -F '\t' '{print $1":"$2,$1,$2,$3,$4,$5,"-",$1":"$2}' OFS='\t' > $outfile.2.tmp
tail -n +2 $outfile | grep -v ^- | awk '$3~/hg18/' | sed 's/hg18//' | awk -F '\t' '{print $1":"$2,$1,$2,$3,$4,$5,"-","-"}' OFS='\t' >> $outfile.2.tmp
sort $outfile.2.tmp | uniq | sort -s -k1,1 > $outfile.2.tmp2
mv $outfile.2.tmp2 $outfile.2.tmp

#run UCSC liftOver for positional map (hg18->hg19)
tail -n +2 $outfile | cut -f1-3 | grep -v ^- | grep hg18 | sed 's/hg18//' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}.hg18.txt
if [ -s ${outfile%.*}.hg18.txt ]; then
    $liftover ${outfile%.*}.hg18.txt $chainfile_hg1819 ${outfile%.*}_liftover_hg18_out_hg19.bed ${outfile%.*}_liftover_unlifted_hg18.bed
    cat ${outfile%.*}_liftover_hg18_out_hg19.bed | awk '{print $4,$2}' OFS='\t' >  ${outfile%.*}_liftover_hg18_to_hg19_map.txt
    cat ${outfile%.*}_liftover_unlifted_hg18.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg18_snps.txt

    #delete unlifted hg18 snps 
    if [ -s ${outfile%.*}_liftover_hg18_unlifted_snps.txt ]; then
        awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' ${outfile%.*}_liftover_hg18_unlifted_snps.txt $outfile.2.tmp > $outfile.2.tmp2
        mv $outfile.2.tmp2 $outfile.2.tmp
    fi

    #update snps with lifted hg19 positions (OUTPUT: CHR:POS CHR POS SNP STUDY TRAIT Hg19 Hg38)
    if [ -s ${outfile%.*}_liftover_hg18_to_hg19_map.txt ]; then
        sort ${outfile%.*}_liftover_hg18_to_hg19_map.txt | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg19_map.tmp
        join -1 1 -2 1 -a 1 -t '	' $outfile.2.tmp ${outfile%.*}_liftover_hg19_map.tmp | awk -F '\t' '{if(NF==9) print "hg19"$1,$2,$9,$4,$5,$6,$2":"$9,$8; else print $0}' OFS='\t' > $outfile.2.tmp2
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
    cat ${outfile%.*}_liftover_hg19_out_hg38.bed | awk '{print $4,$2}' OFS='\t' >  ${outfile%.*}_liftover_hg19_to_hg38_map.txt
    cat ${outfile%.*}_liftover_unlifted_hg19.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg19_snps.txt

    #delete unlifted hg19 snps 
    if [ -s ${outfile%.*}_liftover_unlifted_hg19_snps.txt ]; then
        rm -f ${outfile%.*}_liftover_unlifted_hg19_snps_1.txt
        rm -f ${outfile%.*}_liftover_unlifted_hg19_be_hg19_snps.txt
        cat ${outfile%.*}_liftover_unlifted_hg19_snps.txt | while read snp; do
            chr=${snp%:*}
            rs=`grep "^$snp	" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
            if [ "$rs" = "" ]; then
                rs=`grep "	$snp$" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
            fi
            if [ "$rs" = "" ]; then
                echo $snp >> ${outfile%.*}_liftover_unlifted_hg19_snps_1.txt
            else
                echo $snp"	"$rs >> ${outfile%.*}_liftover_unlifted_hg19_be_hg19_snps.txt
            fi
        done
        if [ -s  ${outfile%.*}_liftover_unlifted_hg19_snps_1.txt ]; then
            awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' ${outfile%.*}_liftover_unlifted_hg19_snps_1.txt $outfile.2.tmp2 > $outfile.2.tmp3
            mv $outfile.2.tmp3 $outfile.2.tmp2
        fi
    fi
    
    #update snps with lifted hg38 positions (OUTPUT: CHR:POS CHR POS SNP STUDY TRAIT Hg19 Hg38)
    if [ -s ${outfile%.*}_liftover_hg19_to_hg38_map.txt ]; then
        sort ${outfile%.*}_liftover_hg19_to_hg38_map.txt | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg38_map.tmp
        join -1 1 -2 1 -a 1 -t '	' $outfile.2.tmp2 ${outfile%.*}_liftover_hg38_map.tmp | awk -F '\t' '{if(NF==9) print $1,$2,$3,$4,$5,$6,$7,$2":"$9; else print $0}' OFS='\t' | sort -s -k1,1 > $outfile.2.tmp3
        rm -f ${outfile%.*}_liftover_hg38_map.tmp
    fi
fi

n_snps=`cat $outfile.2.tmp2 | wc -l`
n_hg19=`awk -F '\t' '$7==$1' $outfile.2.tmp2 | wc -l`
n_hg38=`awk -F '\t' '$8==$1' $outfile.2.tmp2 | wc -l`
n_hg1938_shared=`awk -F '\t' '$8==$7 {print $8}' $outfile.2.tmp3 | wc -l`
n_hg1938_shared_uniq=`sed 's/:/ /' ${outfile%.*}_liftover_hg19_to_hg38_map.txt | awk '$2==$3' | sort | uniq | wc -l`
n_hg19_lack=`awk -F '\t' '$7=="-" {print $8}' $outfile.2.tmp3 | wc -l`
n_hg19_uniq_lack=`awk -F '\t' '$7=="-" {print $8}' $outfile.2.tmp3 | sort | uniq | wc -l`
n_hg38_lack=`awk -F '\t' '$8=="-" {print $7}' $outfile.2.tmp3 | wc -l`
n_hg38_uniq_lack=`awk -F '\t' '$8=="-" {print $7}' $outfile.2.tmp3 | sort | uniq | wc -l`
n_snps_remain=`cut -f4 $outfile.2.tmp3 | wc -l`
echo "n_snps:$n_snps n_hg19:$n_hg19 n_hg38:$n_hg38"
echo "n_hg19_hg39_shared_pos:$n_hg1938_shared(uniq:$n_hg1938_shared_uniq)"
echo "n_hg19_lack:$n_hg19_lack(uniq:$n_hg19_uniq_lack) n_hg38_lack:$n_hg38_lack(uniq:$n_hg38_uniq_lack)"
echo "n_snps_lack_with_hg38_chr:pos:"`awk -F '\t' '$8=="-" {print "hg19:"$7}' $outfile.2.tmp3 | sort | uniq`
echo "#-----------------------------------------------------------------------#"

#run UCSC liftOver for positional map (hg38->hg19)
awk -F '\t' '$7=="-"' $outfile.2.tmp3 | awk '{print "chr"$2,$3,($3+1),$1}' OFS='\t' > ${outfile%.*}.hg38.txt

if [ -s ${outfile%.*}.hg38.txt ]; then
    $liftover ${outfile%.*}.hg38.txt $chainfile_hg3819 ${outfile%.*}_liftover_hg38_out_hg19.bed ${outfile%.*}_liftover_unlifted_hg38.bed
    cat ${outfile%.*}_liftover_hg38_out_hg19.bed | awk '{print $4,$2}' OFS='\t' > ${outfile%.*}_liftover_hg38_to_hg19_map.txt
    cat ${outfile%.*}_liftover_unlifted_hg38.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg38_snps.txt

    if [ -s ${outfile%.*}_liftover_unlifted_hg38_snps.txt ]; then
        rm -f ${outfile%.*}_liftover_unlifted_hg38_snps_1.txt
        rm -f ${outfile%.*}_liftover_unlifted_hg38_be_hg19_snps.txt
        rm -f ${outfile%.*}_liftover_unlifted_hg38_be_hg38_snps.txt
        cat ${outfile%.*}_liftover_unlifted_hg38_snps.txt | while read snp; do
            echo $snp
            chr=${snp%:*}
            #search hg38
            rs=`grep "^$snp	" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f5`
            if [ "$rs" = "" ]; then
                rs=`grep "	$snp$" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f5`
            fi
            if [ "$rs" = "" ]; then
                rs=`grep "^$snp	" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
                if [ "$rs" = "" ]; then
                    rs=`grep "	$snp$" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
                fi
                #search hg19
                if [ "$rs" != "" ]; then
                    echo $snp"	"$rs >>  ${outfile%.*}_liftover_unlifted_hg38_be_hg19_snps.txt
                else
                    echo $snp >> ${outfile%.*}_liftover_unlifted_hg38_snps_1.txt
                fi
            else
                echo $snp"	"$rs >> ${outfile%.*}_liftover_unlifted_hg38_be_hg38_snps.txt
            fi
        done

        #update snps with lifted hg19 positions (OUTPUT: CHR:POS CHR POS SNP STUDY TRAIT Hg19 Hg38)
        if [ -s ${outfile%.*}_liftover_hg38_to_hg19_map.txt ]; then
            sort ${outfile%.*}_liftover_hg38_to_hg19_map.txt | uniq | sort -s -k1,1 > ${outfile%.*}_liftover_hg19_map.tmp
            awk -F '\t' '$7=="-"' $outfile.2.tmp3 > $outfile.2.tmp31
            join -1 1 -2 1 -a 1 -t '	' $outfile.2.tmp31 ${outfile%.*}_liftover_hg19_map.tmp | awk -F '\t' '{if(NF==9) print $1,$2,$9,$4,$5,$6,$2":"$9,$8; else print $0}' OFS='\t' > $outfile.2.tmp4
            awk -F '\t' '$7!="-"' $outfile.2.tmp3 >> $outfile.2.tmp4
            rm -f $outfile.2.tmp31
            rm -f ${outfile%.*}_liftover_hg19_map.tmp
        else
            cp $outfile.2.tmp3 $outfile.2.tmp4
        fi

        n_hg38to19=$((`awk -F '\t' '$7=="-"' $outfile.2.tmp3 | wc -l`-`awk -F '\t' '$7=="-"' $outfile.2.tmp4 | wc -l`))
        n_hg38to19_map=`cat ${outfile%.*}_liftover_hg38_to_hg19_map.txt | wc -l`
        n_hg38to19_map_uniq=`cat ${outfile%.*}_liftover_hg38_to_hg19_map.txt | sort | uniq | wc -l`
        echo "n_hg38to19:$n_hg38to19(map:$n_hg38to19_map uniq:$n_hg38to19_map_uniq)"
        cp $outfile.2.tmp4 $outfile.2.tmp4.bak

        rm -f ${outfile%.*}.hg19_2.txt
        #update unlifted hg38 snps be hg19
        cat ${outfile%.*}_liftover_unlifted_hg38_be_hg19_snps.txt | while read line; do
            snp=`echo "$line" | cut -f1`
            echo $snp | sed 's/:/ /' | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' >> ${outfile%.*}.hg19_2.txt
            awk -v snp=$snp -F '\t' '{if($7=="-" && $1==snp) print $1,$2,$3,$4,$5,$6,snp,"-"; else print $0}' OFS='\t' $outfile.2.tmp4 > $outfile.2.tmp5
            mv $outfile.2.tmp5 $outfile.2.tmp4
        done
        n_unlift_hg38be19=`diff $outfile.2.tmp4.bak $outfile.2.tmp4 | grep \> | cut -d' ' -f2- | cut -f7 | wc -l`
        n_unlift_hg38be19_uniq=`diff $outfile.2.tmp4.bak $outfile.2.tmp4 | grep \> | cut -d' ' -f2- | cut -f7 | sort | uniq | wc -l`
        n_unlift_hg38be19_map=`cat ${outfile%.*}_liftover_unlifted_hg38_be_hg19_snps.txt | wc -l`
        echo "n_unlift_hg38be19:$n_unlift_hg38be19(uniq:$n_unlift_hg38be19_uniq map:$n_unlift_hg38be19_map)"
        cp $outfile.2.tmp4 $outfile.2.tmp4.bak2

        if [ -s ${outfile%.*}_liftover_unlifted_hg38_snps_1.txt ]; then 
            #try to liftover from hg18 to hg19
            sed 's/:/ /' ${outfile%.*}_liftover_unlifted_hg38_snps_1.txt | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' > ${outfile%.*}.hg18_2.txt
            $liftover ${outfile%.*}.hg18_2.txt $chainfile_hg1819 ${outfile%.*}_liftover_hg18_out_hg19_2.bed ${outfile%.*}_liftover_unlifted_hg18_2.bed
            cat ${outfile%.*}_liftover_hg18_out_hg19_2.bed | awk '{print $4,$2}' OFS='\t' > ${outfile%.*}_liftover_hg18_to_hg19_map_2.txt
            cat ${outfile%.*}_liftover_unlifted_hg18_2.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg38_snps_2.txt
            if [ -s ${outfile%.*}_liftover_unlifted_hg38_snps_2.txt ]; then
                sed 's/:/ /' ${outfile%.*}_liftover_unlifted_hg38_snps_2.txt | awk '{print "chr"$1,$2,($2+1),$1":"$2}' OFS='\t' >> ${outfile%.*}.hg19_2.txt
            fi

            #update lifted hg18 snp postions to be hg19
            cat ${outfile%.*}_liftover_hg18_to_hg19_map_2.txt | while read line; do
                snp=`echo "$line" | cut -f1`
                pos=`echo "$line" | cut -f2`
                snp2=${snp%:*}":"$pos
                awk -v snp=$snp -v snp2=$snp2 -v pos=$pos -F '\t' '{if($7=="-" && $1==snp) print $1,$2,pos,$4,$5,$6,snp2,"-"; else print $0}' OFS='\t' $outfile.2.tmp4 > $outfile.2.tmp5
                echo "chr"${snp%:*}"	"$pos"	"$((pos+1))"	"$snp >> ${outfile%.*}.hg19_2.txt
                mv $outfile.2.tmp5 $outfile.2.tmp4
            done

            n_hg18to19_2=`diff $outfile.2.tmp4.bak2 $outfile.2.tmp4 | grep \> | cut -d' ' -f2- | cut -f7 | sort | uniq | wc -l`
            n_hg18to19_map_2=`cat ${outfile%.*}_liftover_hg18_to_hg19_map_2.txt | wc -l`
            n_hg18to19_map_uniq_2=`cat ${outfile%.*}_liftover_hg18_to_hg19_map_2.txt | sort | uniq | wc -l`
            echo "n_hg18to19_2:$n_hg18to19_2(map:$n_hg18to19_map_2 uniq:$n_hg18to19_map_uniq_2)"
            cp $outfile.2.tmp4 $outfile.2.tmp4.bak3

            #try to liftover updated snps from hg19 to hg38
            $liftover ${outfile%.*}.hg19_2.txt $chainfile_hg1938 ${outfile%.*}_liftover_hg19_out_hg38_2.bed ${outfile%.*}_liftover_unlifted_hg19_2.bed
            cat ${outfile%.*}_liftover_hg19_out_hg38_2.bed | awk '{print $4,$2}' OFS='\t' > ${outfile%.*}_liftover_hg19_to_hg38_map_2.txt
            cat ${outfile%.*}_liftover_unlifted_hg19_2.bed | grep ^chr | awk '{print $4}' | sort -n | uniq > ${outfile%.*}_liftover_unlifted_hg38_snps_3.txt

            #update lifted hg19 snp postions to be hg38
            cat ${outfile%.*}_liftover_hg19_to_hg38_map_2.txt | while read line; do
                snp=`echo "$line" | cut -f1`
                pos=`echo "$line" | cut -f2`
                snp2=${snp%:*}":"$pos
                awk -v snp=$snp -v snp2=$snp2 -v pos=$pos -F '\t' '{if($1==snp && $7!="-") print $1,$2,$3,$4,$5,$6,$7,snp2; else if($1==snp && $7=="-") print $1,$2,$3,$4,$5,$6,snp,snp2; else print $0}' OFS='\t' $outfile.2.tmp4 > $outfile.2.tmp5
                mv $outfile.2.tmp5 $outfile.2.tmp4
            done
            n_hg19to38_2=`diff $outfile.2.tmp4.bak3 $outfile.2.tmp4 | grep \> | cut -d' ' -f2- | cut -f8 | wc -l`
            n_hg19to38_uniq_2=`diff $outfile.2.tmp4.bak3 $outfile.2.tmp4 | grep \> | cut -d' ' -f2- | cut -f8 | sort | uniq | wc -l`
            n_hg19to38_map_2=`cat ${outfile%.*}_liftover_hg19_to_hg38_map_2.txt | wc -l`
            n_hg19to38_map_uniq_2=`cat ${outfile%.*}_liftover_hg19_to_hg38_map_2.txt | sort | uniq | wc -l`
            echo "n_hg19to38_2:$n_hg19to38_2(uniq:$n_hg19to38_uniq_2 map:$n_hg19to38_map_2 uniq:$n_hg19to38_map_uniq_2)"
            cp $outfile.2.tmp4 $outfile.2.tmp4.bak4
        fi

        #delete unlifted hg38 snps 
        if [ -s  ${outfile%.*}_liftover_unlifted_hg38_snps_3.txt ]; then
            awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' ${outfile%.*}_liftover_unlifted_hg38_snps_3.txt $outfile.2.tmp4 > $outfile.2.tmp5
            mv $outfile.2.tmp5 $outfile.2.tmp4
        fi
    fi

    awk -F '\t' '{if($7!="-"&&$2":"$3==$7) print $1,$2,$3,$4,$5,$6,$7,$8,"hg19"; else if($8!="-"&&$2":"$3==$8) print $1,$2,$3,$4,$5,$6,$7,$8,"hg38"}' OFS='\t' $outfile.2.tmp4 > $outfile.2.tmp5
fi
awk -F '\t' '$4!~"rs"' $outfile.2.tmp5 | while read line; do
    snp=`echo "$line" | awk -F '\t' '{print $2":"$3}'`
    chr=`echo "$line" | cut -f2`
    tag=`echo "$line" | cut -f9`
    if [ "$tag" = "hg19" ]; then
        rs=`grep "^$snp	" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
        if [ "$rs" = "" ]; then
            rs=`grep "	$snp$" $(dirname $outfile)/hg19/hg19_snps.chr$chr.csv | cut -f5`
        fi
    elif [ "$tag" = "hg38" ]; then
        rs=`grep "^$snp	" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f5`
        if [ "$rs" = "" ]; then
            rs=`grep "	$snp$" $(dirname $outfile)/hg38/hg38_snps.chr$chr.csv | cut -f5`
        fi
    fi
    if [ "$rs" != "" ]; then
        echo "$line" | awk -v rs=$rs -F '\t' '{print $1,$2,$3,rs,$5,$6,$7,$8,$9}'
    else
        echo "$line"
    fi
    break
done
fi
#rm -f $outfile*.tmp*
#-------------------------------------------------------------------------#
