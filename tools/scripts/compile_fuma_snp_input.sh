#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to compile fuma input snp list across multiple pleioFDR
# analyses.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh compile_fuma_snp_input.sh compile_config outfolder"
  echo "Arguments: compile_config - compile config file"
  echo "           outfolder - output folder"
  echo "Example:   sh compile_fuma_snp_input.sh ./mood_psych_compile_config.txt outfolder"
  exit 0
fi
#------------------------------------------------------------------------#
compile_config=$1
outfolder=$2

if [ ! -f $compile_config ]; then
    echo "NOTE: $compile_config doesn't exist" 
    exit 1
fi

source $compile_config

if [ -z "$taglist" ]; then
    echo "NOTE: please set 'taglist' variable"
    exit 1
fi

if [ -z "$sumstatlist" ]; then
    echo "NOTE: please set 'sumstatlist' variable"
    exit 1
fi

if [ -z $pleiofdr_outfolder ]; then
    echo "NOTE: please set 'pleiofdr_outfolder' variable"
    exit 1
fi

if [ -z "$complex_regions" ]; then
    echo "NOTE: please set 'complex_regions' variable"
    exit 1
fi

if [ -z $rev_flag ]; then
    echo "NOTE: please set 'rev_flag' variable"
    exit 1
fi

if [ ! -d $pleiofdr_outfolder ]; then
    echo "NOTE: $pleiofdr_outfolder doesn't exist"
    exit 1
fi
#------------------------------------------------------------------------#

tagA=`echo $taglist | awk '{print $1}'`
tagB=`echo $taglist | awk '{print $2}' | awk -F ':' '{print $1}'`
IFS='#' read -r -a sumstat_array <<< "$sumstatlist"

n_tag=`echo $taglist | sed 's/:/ /' | sed 's/,/ /g' | tr -s ' ' | sed 's/ /\n/g' | wc -l`
n_sumstat=${#sumstat_array[@]}
if [ $((n_tag-1)) -ne $n_sumstat ]; then 
   echo "NOTE: the number of tags doesn't match the number of sumstats"
   exit 1
fi

if [ "$rev_flag" = "Y" ]; then
    tag11=$tagB
    tag22=$tagA
else
    tag11=$tagA
    tag22=$tagB
fi

rm -f $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
n=1
for tag in `echo $taglist | awk -F ':' '{print $2}' | sed 's/,/\n/g'`; do
    if [ "$rev_flag" = "Y" ]; then
        tag1=$tag
        tag2=$tagA
        sm1=`echo ${sumstat_array[$n]}`
        sm2=`echo ${sumstat_array[0]}`
    else
        tag1=$tagA
        tag2=$tag
        sm1=`echo ${sumstat_array[0]}`
        sm2=`echo ${sumstat_array[$n]}`
    fi

    n=$((n+1))

    echo $tag1 $tag2
    echo $sm1 $sm2

    conjfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_conjfdr

    if [ ! -d $conjfdr ]; then
        echo "NOTE: $conjfdr doesn't exist, try to switch rev_flag"
        continue
    fi

    condfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_condfdr
    if [ ! -d $condfdr ]; then
        echo "NOTE: $condfdr doesn't exist"
        continue
    fi

    condfdr2=$pleiofdr_outfolder/${sm2}_vs_${sm1}_condfdr
    if [ ! -d $condfdr2 ]; then
        echo "NOTE: $condfdr2 doesn't exist"
        continue
    fi

    sh $tab/fdr_fuma_snp_input.sh $conjfdr/conj.result.clump.snps.csv 0.1 0.6 $outfolder/${tag1}_vs_${tag2}_conj_snps.txt "$complex_regions"
    sh $tab/fdr_fuma_snp_input.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $outfolder/${tag1}_vs_${tag2}_cond_snps.txt "$complex_regions"
    sh $tab/fdr_fuma_snp_input.sh $condfdr2/cond.result.clump.snps.csv 0.1 0.6 $outfolder/${tag2}_vs_${tag1}_cond_snps.txt "$complex_regions"
    cat $outfolder/${tag1}_vs_${tag2}_cond_snps.txt $outfolder/${tag2}_vs_${tag1}_cond_snps.txt | grep -v snp >> $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
    rm -f $outfolder/${tag1}_vs_${tag2}_cond_snps.txt $outfolder/${tag2}_vs_${tag1}_cond_snps.txt
done

if [ -s $outfolder/${tag11}_vs_${tag22}_cond_snps.txt ]; then
    echo 'snp	pval' > $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp
    cat $outfolder/${tag11}_vs_${tag22}_cond_snps.txt | cut -f1 | sort -s | uniq |  awk '{print $1,1e-10}' OFS='\t' >>  $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp
    mv $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
fi
