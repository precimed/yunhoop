account=nn9114k
dbsnp_merge_map_file=$1
target_rs_list=$2
target_merge_map_file=$3
chunksize=$4

awk -F '\t' 'NR==FNR {D[$1]++; next} ($1 in D)' $target_rs_list $dbsnp_merge_map_file | cut -f1-2 > $target_merge_map_file
awk -F '\t' 'NR==FNR {D[$1]++; next} !($2 in D)' $dbsnp_merge_map_file $target_merge_map_file | cut -f1-2 > ${target_merge_map_file%.*}_B_0.txt
if [ `cat ${target_merge_map_file%.*}_B_0.txt | wc -l` -eq `cat $target_merge_map_file | wc -l` ]; then
    rm -f ${target_merge_map_file%.*}_B_0.txt
    echo "ready"
    exit 0
fi
awk -F '\t' 'NR==FNR {D[$1]++; next} ($2 in D)' $dbsnp_merge_map_file $target_merge_map_file | cut -f1-2 > ${target_merge_map_file%.*}_A.txt
checksize=`cat ${target_merge_map_file%.*}_A.txt | wc -l`
if [ -z $chunksize ]; then
    chunksize=1000
fi
chunknum=`echo $checksize $chunksize | awk '{print $1/$2}' | awk '{print $1=$1==int($1)?int($1):int($1)+1}'`
if [ $chunknum -lt 10 ]; then
    len=1
elif [ $chunknum -ge 10 ] && [ $chunknum -le 99 ]; then
    len=2
elif [ $chunknum -ge 100 ] && [ $chunknum -le 999 ]; then
    len=3
elif [ $chunknum -ge 1000 ] && [ $chunknum -le 9999 ]; then
    len=4
else
    len=5
fi
echo $checksize $chunksize $len

split -l $chunksize --numeric-suffixes=1 --suffix-length=$len --additional-suffix=".txt" ${target_merge_map_file%.*}_A.txt ${target_merge_map_file%.*}_A_
for i in ${target_merge_map_file%.*}_A_*.txt; do
    suffix=${i%.*}
    suffix=${suffix##*_}
    echo "#!/bin/sh" > $(dirname $target_merge_map_file)/build_snp_merge_map_$suffix.sh
    echo "sh $(dirname $0)/track_snp_merge.sh $dbsnp_merge_map_file $i ${target_merge_map_file%.*}_B_$suffix.txt &> ${target_merge_map_file%.*}_A_$suffix.log" >> $(dirname $target_merge_map_file)/build_snp_merge_map_$suffix.sh
    chmod +x $(dirname $target_merge_map_file)/build_snp_merge_map_$suffix.sh
    if [ $checksize -le 200 ]; then
        sh $(dirname $target_merge_map_file)/build_snp_merge_map_$suffix.sh
    else
        srun -A $account -n1 --mem 4G -t 12:00:00 $(dirname $target_merge_map_file)/build_snp_merge_map_$suffix.sh &
    fi
done
wait

cat ${target_merge_map_file%.*}_B_*.txt | sort -s -k1,1 > ${target_merge_map_file%.*}_B.txt
awk -F '\t' '{print $2,$1}' OFS='\t' ${target_merge_map_file%.*}_B.txt | uniq -c -f1 | awk '{print $3,$1}' OFS='\t' > ${target_merge_map_file%.*}_Bn.txt
join -1 1 -2 1 -t '	' ${target_merge_map_file%.*}_B.txt ${target_merge_map_file%.*}_Bn.txt | awk -F '\t' '{print $1,$2,$3}' OFS='\t' > ${target_merge_map_file%.*}_2.txt

if [ `cat ${target_merge_map_file%.*}_A*.log | wc -l` -eq 0 ]; then
    rm -f ${target_merge_map_file%.*}_A*.txt ${target_merge_map_file%.*}_B*.txt
    rm -f $(dirname $target_merge_map_file)/build_snp_merge_map_*.sh
    rm -f ${target_merge_map_file%.*}_A*.log
    mv ${target_merge_map_file%.*}_2.txt ${target_merge_map_file%.*}_full.txt
fi

#only take first entry of each merging group
#cut -f1 ${target_merge_map_file%.*}_full.txt | uniq | while read snp; do awk -v snp=$snp '$1==snp {print; exit}' ${target_merge_map_file%.*}.txt; done
