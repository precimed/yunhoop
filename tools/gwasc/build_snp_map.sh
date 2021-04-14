account=nn9114k
dbsnp_folder=$1
target_rs_list=$2
target_map_folder=$3

if [ ! -d $dbsnp_folder/json ]; then
    echo "please place your downloaded json files in $dbsnp_folder/json folder"
    exit 1
fi
if [ ! -d $dbsnp_folder/vcf ]; then
    echo "please place your downloaded vcf files in $dbsnp_folder/vcf folder"
    exit 1
fi

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M; do
    continue
    echo "#!/bin/sh" > $target_map_folder/build_snp_map_A_chr$i.sh
    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/json/dbsnp_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_A_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)" >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo 'if [ `grep GRCh38 '"$target_map_folder/gwasc_A_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) | wc -l"'` -gt 0 ]; then' >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo "    grep GRCh38 $target_map_folder/gwasc_A_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) > $target_map_folder/gwasc_A_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)" >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo "fi" >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo 'if [ `grep GRCh37 '"$target_map_folder/gwasc_A_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) | wc -l"'` -gt 0 ]; then' >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo "    grep GRCh37 $target_map_folder/gwasc_A_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) > $target_map_folder/gwasc_A_b37_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)" >> $target_map_folder/build_snp_map_A_chr$i.sh
    echo "fi" >> $target_map_folder/build_snp_map_A_chr$i.sh
    chmod +x $target_map_folder/build_snp_map_A_chr$i.sh
    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_A_chr$i.sh &
done

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M T W; do
    continue
    echo "#!/bin/sh" > $target_map_folder/build_snp_map_B_b37_chr$i.sh
    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/vcf/dbsnp_b37_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_B_b37_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)" >> $target_map_folder/build_snp_map_B_b37_chr$i.sh
    chmod +x $target_map_folder/build_snp_map_B_b37_chr$i.sh
    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_B_b37_chr$i.sh & 

    echo "#!/bin/sh" > $target_map_folder/build_snp_map_B_b38_chr$i.sh
    echo "awk -F '\t' 'NR==FNR "'{D[$1]++; next} ($1 in D)'"' $target_rs_list $dbsnp_folder/vcf/dbsnp_b38_chr$i.txt | sort -s -k1,1 > $target_map_folder/gwasc_B_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)" >> $target_map_folder/build_snp_map_B_b38_chr$i.sh
    chmod +x $target_map_folder/build_snp_map_B_b38_chr$i.sh
    srun -A $account -n1 --mem 4G -t 1:00:00 $target_map_folder/build_snp_map_B_b38_chr$i.sh & 
done
wait

rm -f $target_map_folder/build_snp_map_A_chr*.sh
rm -f $target_map_folder/build_snp_map_B_b*_chr*.sh

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    join -1 1 -2 1 -t'	' $target_map_folder/gwasc_A_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) $target_map_folder/gwasc_B_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) | awk -F '\t' '$2!=$10' > $target_map_folder/gwasc_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt)
    join -1 1 -2 1 -t'	' $target_map_folder/gwasc_A_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) $target_map_folder/gwasc_B_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) | awk -F '\t' '$2!=$10' | cut -f1-2,9-10 > $target_map_folder/gwasc_b38_$(basename $dbsnp_folder/json/dbsnp_chr${i}_diff.txt)
    join -1 1 -2 1 -t'	' $target_map_folder/gwasc_A_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) $target_map_folder/gwasc_B_b38_$(basename $dbsnp_folder/json/dbsnp_chr$i.txt) | awk -F '\t' '$2!=$10 && $9!=$10' | cut -f1-2,9-10 > $target_map_folder/gwasc_b38_$(basename $dbsnp_folder/json/dbsnp_chr${i}_diff2.txt)
done
