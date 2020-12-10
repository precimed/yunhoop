#!/bin/bash
#--------------------------- Description ---------------------------------#
# This script converts a general mostest csv to sumstat, but WITHOUT ANY
# WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 5 ]; then
    echo "Usage:     sh sumstat.sh python_convert misc rawfile outname outfolder"
    echo "Arguments: python_convert - home folder to python_convert
                     misc - folder containing reference files of 9545380.ref and ega.grch37.all.haplotypes.m.v2.ref
                     rawfile - file of raw data
                     outname - output name
                     outfolder - output folder"
    exit 0
fi
#-------------------------------------------------------------------------#

module purge
module load Python/3.7.4-GCCcore-8.3.0
set -o errexit

python_convert=$1
misc=$2
raw=$3
name=$4
sumstat=$5

mkdir -p $sumstat/STD
mkdir -p $sumstat/TMP
mkdir -p $sumstat/TMP/csv
mkdir -p $sumstat/TMP/zscore
mkdir -p $sumstat/TMP/qc
mkdir -p $sumstat/TMP/lift
mkdir -p $sumstat/TMP/variantid
mkdir -p $sumstat/TMP/mat_9545380
mkdir -p $sumstat/TMP/neff
mkdir -p $sumstat/TMP/nomhc

python $python_convert/sumstats.py csv \
    --force  \
    --out $sumstat/TMP/csv/$name.sumstats \
    --log $sumstat/STD/$name.log \
    --auto  \
    --sumstats $raw

gzip $sumstat/TMP/csv/$name.sumstats

python $python_convert/sumstats.py zscore \
    --force  \
    --sumstats $sumstat/TMP/csv/$name.sumstats.gz \
    --log-append  \
    --out $sumstat/TMP/zscore/$name.sumstats \
    --log $sumstat/STD/$name.log

gzip $sumstat/TMP/zscore/$name.sumstats

python $python_convert/sumstats.py qc \
    --force  \
    --log $sumstat/STD/$name.log \
    --sumstats $sumstat/TMP/zscore/$name.sumstats.gz \
    --max-or 1e+37 \
    --log-append  \
    --fix-dtype-cols ['BP', 'CHR', 'PVAL'] \
    --dropna-cols ['A1', 'A2'] \
    --out $sumstat/TMP/qc/$name.sumstats

gzip $sumstat/TMP/qc/$name.sumstats
cp $sumstat/TMP/qc/$name.sumstats.gz $sumstat/TMP/lift/$name.sumstats.gz

python $python_convert/sumstats.py variantid \
    --force  \
    --log $sumstat/STD/$name.log \
    --sumstats $sumstat/TMP/lift/$name.sumstats.gz \
    --out $sumstat/TMP/variantid/$name.sumstats \
    --ref $misc/ega.grch37.all.haplotypes.m.v2.ref \
    --log-append

gzip $sumstat/TMP/variantid/$name.sumstats
cp $sumstat/TMP/variantid/$name.sumstats.gz $sumstat/STD/$name.sumstats.gz

python $python_convert/sumstats.py mat \
    --force  \
    --sumstats $sumstat/STD/$name.sumstats.gz \
    --out $sumstat/TMP/mat_9545380/$name.mat \
    --keep-all-cols  \
    --ref $misc/9545380_ref/9545380.ref

python $python_convert/sumstats.py neff \
    --force  \
    --sumstats $sumstat/STD/$name.sumstats.gz \
    --drop  \
    --out $sumstat/TMP/neff/$name.sumstats

gzip $sumstat/TMP/neff/$name.sumstats

python $python_convert/sumstats.py qc \
    --force  \
    --sumstats $sumstat/TMP/neff/$name.sumstats.gz \
    --exclude-ranges 6:26000000-34000000 \
    --out $sumstat/TMP/nomhc/$name.sumstats

gzip $sumstat/TMP/nomhc/$name.sumstats
