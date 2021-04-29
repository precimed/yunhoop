if [ $# -lt 1 ]; then
    echo "sh rsjson.sh chr_rs_list[chr1:rs1,chr2:rs2,chr3:rs3]"
    exit 0
fi

chr_rs_list=$1
account=nn9114k

for ln in `echo $chr_rs_list | sed 's/,/ /g'`; do
    chr=${ln%:*}
    rs=${ln#*:}
    rs=${rs/rs/}
    echo $chr $rs
    echo "#!/bin/sh" > tmp_rs$rs.sh
    echo "sh $db/readjson_rs.sh refsnp-chr$chr.json.bz2 $rs > ${chr}_rs$rs.txt" >> tmp_rs$rs.sh
    chmod +x tmp_rs$rs.sh
    srun -A $account -n1 --mem 4G -t 2:00:00 tmp_rs$rs.sh &
done
wait
rm -f tmp_rs*.sh
