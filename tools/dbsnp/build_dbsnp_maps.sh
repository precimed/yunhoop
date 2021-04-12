input_folder=$1
out_map_folder=$2
account=nn9114k
host=https://ftp.ncbi.nih.gov/snp/redesign/archive/b153/JSON

#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
#    wget --continue --tries=0 $host/refsnp-chr$i.json.bz2 -O $input_folder/refsnp-chr$i.json.bz2
#    wget --continue --tries=0 $host/refsnp-chr$i.json.bz2.md5 -O $input_folder/refsnp-chr$i.json.bz2.md5
#done
#for i in merged withdrawn unsupported nosnppos other; do
#    wget --continue --tries=0 $host/refsnp-$i.json.bz2 -O $input_folder/refsnp-$i.json.bz2
#    wget --continue --tries=0 $host/refsnp-$i.json.bz2.md5 -O $input_folder/refsnp-$i.json.bz2.md5
#    md5sum -c $input_folder/refsnp-$i.json.bz2.md5
#done

#------------------------build misc maps------------------------#
python rsjson_misc_map.py -m refsnp-merged.json -n refsnp-nosnppos.json -u refsnp-unsupported.json -w refsnp-withdrawn.json -o refsnp-other.json -t $out_map_folder &> $out_map_folder/refsnp_misc.log &
#---------------------------------------------------------------#


#------------------------build chr maps------------------------#
#for i in 5 1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 X Y MT; do
#    md5sum -c refsnp-chr$i.json.bz2.md5
#done
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 X Y MT; do
#for i in 11 12 13 14; do
#for i in 15 16 17 18; do
#for i in 19 20 21 22; do
#for i in X Y MT; do
for i in MT; do
    chr=$i
    if [ $i = "MT" ]; then
        chr='M2'
    fi
    python rsjson_chr_map.py -i $input_folder/refsnp-chr$i.json.bz2 -c $chr -o $out_map_folder/dbsnp_chr$chr.txt &> $out_map_folder/dbsnp_chr$chr.log &
done
#---------------------------------------------------------------#
