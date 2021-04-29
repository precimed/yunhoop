#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to 

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 3 ]; then
  echo "Usage: sh build_dbsnp_maps.sh dbsnp_folder account flag_chr flag_misc"
  echo "Arguments: dbsnp_folder - top folder to place dbsnp data"
  echo "           account - slurm account for srun"
  echo "           flag_chr - flag to build chr map [Y/N]"
  echo "           flag_misc - flag to build misc map [Y/N]"
  echo "Example: sh build_dbsnp_maps.sh $lf/ncbi/b153 nn9114k Y N"
  exit 0
fi
#-------------------------------------------------------------------------#

account=nn9114k
dbsnp_folder=$1
flag_chr_A=$2
flag_chr_B=$3
flag_misc=$4

#------------------------build misc maps------------------------#
if [ "$flag_misc" = "Y" ]; then
    echo "...misc"
    echo "#!/bin/sh" > $dbsnp_folder/json/build_dbsnp_misc_map.sh
    echo "python $(dirname $0)/rsjson_misc_map.py -m $dbsnp_folder/json/refsnp-merged.json -n $dbsnp_folder/json/refsnp-nosnppos.json -u $dbsnp_folder/json/refsnp-unsupported.json -w $dbsnp_folder/json/refsnp-withdrawn.json -o $dbsnp_folder/json/refsnp-other.json -t $dbsnp_folder/json &> $dbsnp_folder/json/refsnp_misc.log" >> $dbsnp_folder/json/build_dbsnp_misc_map.sh
    chmod +x $dbsnp_folder/json/build_dbsnp_misc_map.sh
    srun -A $account -n1 --mem 4G -t 5:00:00 $dbsnp_folder/json/build_dbsnp_misc_map.sh &
fi
wait
#---------------------------------------------------------------#

#------------------------build chr maps------------------------#
if [ "$flag_chr_B" = "Y" ]; then
    echo "...vcf"
    for j in 25 38; do
        id=$j
        if [ $j = "25" ]; then
            id='37'
        fi
        echo "#!/bin/sh" > $dbsnp_folder/vcf/build_dbsnp_map_b$id.sh
        echo "zcat $dbsnp_folder/vcf/GCF_000001405.$j.gz | tail -n +38 | cut -f1-5,8 | sed 's/RS=.*VC=//' | sed 's/;.*//' > $dbsnp_folder/vcf/dbsnp_b$id.txt" >> $dbsnp_folder/vcf/build_dbsnp_map_b$id.sh
        chmod +x $dbsnp_folder/vcf/build_dbsnp_map_b$id.sh
        srun -A $account -n1 --mem 4G -t 3:00:00 $dbsnp_folder/vcf/build_dbsnp_map_b$id.sh &
    done
    wait

    for id in 37 38; do
        for ((i=1; i<=27; i++)); do
            chr=$i
            if [ $i = "23" ]; then
                chr="X"
            fi
            if [ $i = "24" ]; then
                chr="Y"
            fi
            if [ "$i" = "25" ]; then
                chr="M"
                seq_id="NC_012920"
            elif [ "$i" = "26" ]; then
                chr="T"
                seq_id="NT_"
            elif [ "$i" = "27" ]; then
                chr="W"
                seq_id="NW_"
            else
                seq_id=$((1000000+i))
                seq_id="NC_"${seq_id:1:6}
            fi
            echo "#!/bin/sh" > $dbsnp_folder/vcf/build_dbsnp_map_b${id}_chr$chr.sh
            echo "grep ^$seq_id $dbsnp_folder/vcf/dbsnp_b$id.txt | awk -v chr=$chr -F '\t' '{print "'$3,chr":"$2,chr,$2,$4,$5,$1,$6}'"' OFS='\t' > $dbsnp_folder/vcf/dbsnp_b${id}_chr$chr.txt" >> $dbsnp_folder/vcf/build_dbsnp_map_b${id}_chr$chr.sh
            chmod +x $dbsnp_folder/vcf/build_dbsnp_map_b${id}_chr$chr.sh
            srun -A $account -n1 --mem 4G -t 3:00:00 $dbsnp_folder/vcf/build_dbsnp_map_b${id}_chr$chr.sh &
        done
    done
fi

if [ "$flag_chr_A" = "Y" ]; then
    echo "...json"
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
        chr=$i
        if [ $i = "MT" ]; then
            chr='M'
        fi
        echo "#!/bin/sh" > $dbsnp_folder/json/build_dbsnp_map_chr$chr.sh
        echo "module load Python/3.7.4-GCCcore-8.3.0" >> $dbsnp_folder/json/build_dbsnp_map_chr$chr.sh
        echo "python $(dirname $0)/rsjson_chr_map.py -i $dbsnp_folder/json/refsnp-chr$i.json.bz2 -c $chr -o $dbsnp_folder/json/dbsnp_chr$chr &> $dbsnp_folder/json/dbsnp_chr$chr.log" >> $dbsnp_folder/json/build_dbsnp_map_chr$chr.sh
        chmod +x $dbsnp_folder/json/build_dbsnp_map_chr$chr.sh
        srun -A $account -n1 --mem 4G -t 6:00:00 $dbsnp_folder/json/build_dbsnp_map_chr$chr.sh &
    done
fi
wait

if [ `cat $dbsnp_folder/json/dbsnp_chr*.log | wc -l` -eq 0 ]; then
    rm -f $dbsnp_folder/json/dbsnp_chr*.log
    rm -f $dbsnp_folder/json/build_dbsnp_map_chr*.sh
fi
rm -f $dbsnp_folder/json/build_dbsnp_misc_map.sh
rm -f $dbsnp_folder/vcf/build_dbsnp_map_b*.sh
rm -f $dbsnp_folder/vcf/build_dbsnp_map_b*_chr*.sh
#---------------------------------------------------------------#
