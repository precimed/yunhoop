account=nn9114k
dbsnp_folder=$1
target_rs_list=$2
target_map_folder=$3

#if [ ! -d $dbsnp_folder/json ]; then
#    echo "please place your downloaded json files in $dbsnp_folder/json folder"
#    exit 1
#fi
if [ ! -d $dbsnp_folder/vcf ]; then
    echo "please place your downloaded vcf files in $dbsnp_folder/vcf folder"
    exit 1
fi

#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
#    echo "#!/bin/sh" > $target_map_folder/build_snp_map_A_chr$i.sh
#    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/json/dbsnp_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_A_dbsnp_chr$i.txt" >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo 'if [ `grep GRCh38 '"$target_map_folder/gwasc_A_dbsnp_chr$i.txt | wc -l"'` -gt 0 ]; then' >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo "    grep GRCh38 $target_map_folder/gwasc_A_dbsnp_chr$i.txt > $target_map_folder/gwasc_A_b38_dbsnp_chr$i.txt" >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo "fi" >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo 'if [ `grep GRCh37 '"$target_map_folder/gwasc_A_dbsnp_chr$i.txt | wc -l"'` -gt 0 ]; then' >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo "    grep GRCh37 $target_map_folder/gwasc_A_dbsnp_chr$i.txt > $target_map_folder/gwasc_A_b37_dbsnp_chr$i.txt" >> $target_map_folder/build_snp_map_A_chr$i.sh
#    echo "fi" >> $target_map_folder/build_snp_map_A_chr$i.sh
#    chmod +x $target_map_folder/build_snp_map_A_chr$i.sh
#    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_A_chr$i.sh &
#done

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M T W; do
    echo "#!/bin/sh" > $target_map_folder/build_snp_map_B_b37_chr$i.sh
    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/vcf/dbsnp_b37_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_B_b37_dbsnp_chr$i.txt" >> $target_map_folder/build_snp_map_B_b37_chr$i.sh
    chmod +x $target_map_folder/build_snp_map_B_b37_chr$i.sh
    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_B_b37_chr$i.sh & 

    echo "#!/bin/sh" > $target_map_folder/build_snp_map_B_b38_chr$i.sh
    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/vcf/dbsnp_b38_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_B_b38_dbsnp_chr$i.txt" >> $target_map_folder/build_snp_map_B_b38_chr$i.sh
    chmod +x $target_map_folder/build_snp_map_B_b38_chr$i.sh
    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_B_b38_chr$i.sh & 
done
wait

#rm -f $target_map_folder/build_snp_map_A_chr*.sh
#rm -f $target_map_folder/build_snp_map_B_b*_chr*.sh
#rm -f $target_map_folder/gwasc_dbsnp_diff_full.txt
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
#    join -1 1 -2 1 -t'	' $target_map_folder/gwasc_A_b38_dbsnp_chr$i.txt $target_map_folder/gwasc_B_b38_dbsnp_chr$i.txt | awk -F '\t' '$2!=$10' > $target_map_folder/gwasc_b38_dbsnp_chr$i.txt
#    join -1 1 -2 1 -t'	' $target_map_folder/gwasc_b38_dbsnp_chr$i.txt $target_map_folder/gwasc_B_b37_dbsnp_chr$i.txt >> $target_map_folder/gwasc_dbsnp_diff_full.txt
#    rm -f $target_map_folder/gwasc_b38_dbsnp_chr$i.txt
#done
#cut -f1-3,9-10,16 $target_map_folder/gwasc_dbsnp_diff_full.txt > $target_map_folder/gwasc_dbsnp_diff.txt
#awk '$4!=$5' $target_map_folder/gwasc_dbsnp_diff.txt > $target_map_folder/gwasc_dbsnp_diff2.txt

#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
#    awk -F '\t' 'NR==FNR {D[$1]++; next} !($1 in D)' $target_map_folder/gwasc_dbsnp_diff.txt $target_map_folder/gwasc_B_b38_dbsnp_chr$i.txt > $target_map_folder/gwasc_C_b38_dbsnp_chr$i.txt
#    awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' $target_map_folder/gwasc_dbsnp_diff.txt $target_map_folder/gwasc_B_b38_dbsnp_chr$i.txt > $target_map_folder/gwasc_B_b38_dbsnp_chr$i.tmp
#    sort -s -k1,1 $target_map_folder/gwasc_dbsnp_diff.txt > $target_map_folder/gwasc_dbsnp_diff.tmp
#    join -1 1 -2 1 -t '	' $target_map_folder/gwasc_B_b38_dbsnp_chr$i.tmp $target_map_folder/gwasc_dbsnp_diff.tmp | awk -F '\t' '{if(NF>7) print $1,$8,$3,$4,$5,$6,$7; else print $1,$2,$3,$4,$5,$6,$7}' OFS='\t' >>  $target_map_folder/gwasc_C_b38_dbsnp_chr$i.txt
#    rm -f $target_map_folder/gwasc_dbsnp_diff.tmp $target_map_folder/gwasc_B_b38_dbsnp_chr$i.tmp
#done
