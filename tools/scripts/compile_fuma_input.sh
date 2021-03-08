#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to compile fuma input snp list across multiple pleioFDR
# analyses.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh compile_fuma_input.sh compile_config outfolder"
  echo "Arguments: compile_config - compile config file"
  echo "           outfolder - output folder"
  echo "Example:   sh compile_fuma_input.sh mood_psych_compile_config.txt outfolder"
  exit 0
fi
#------------------------------------------------------------------------#
compile_config=$1
outfolder=$2

if [ ! -f $compile_config ]; then
    echo "$compile_config doesn't exist" 
    exit 1
fi

source $compile_config

if [ -z "$taglist" ]; then
    echo "please set 'taglist' variable"
    exit 1
fi

if [ -z "$sumstatlist" ]; then
    echo "please set 'sumstatlist' variable"
    exit 1
fi

if [ -z $pleiofdr_outfolder ]; then
    echo "please set 'pleiofdr_outfolder' variable"
    exit 1
fi

if [ -z "$complex_regions" ]; then
    echo "please set 'complex_regions' variable"
    exit 1
fi

if [ -z $rev_flag ]; then
    echo "please set 'rev_flag' variable"
    exit 1
fi

if [ ! -d $pleiofdr_outfolder ]; then
    echo "$pleiofdr_outfolder doesn't exist"
    exit 1
fi
#------------------------------------------------------------------------#

tagA=`echo $taglist | awk '{print $1}'`
tagB=`echo $taglist | awk '{print $2}'`
IFS='#' read -r -a sumstat_array <<< "$sumstatlist"

if [ "$rev_flag" = "Y" ]; then
    tag11=$tagB
    tag22=$tagA
else
    tag11=$tagA
    tag22=$tagB
fi

rm -f $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
n=1
for tag in `echo $taglist | awk '{for(i=3;i<=NF;i++) printf "%s ", $i}'`; do
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
    echo $tag1 $tag2
    echo $sm1 $sm2

    conjfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_conjfdr

    if [ ! -d $conjfdr ]; then
        echo "$conjfdr doesn't exist, try to switch rev_flag"
        continue
    fi

    condfdr=$pleiofdr_outfolder/${sm1}_vs_${sm2}_condfdr
    if [ ! -d $condfdr ]; then
        echo "$condfdr doesn't exist"
        continue
    fi

    condfdr2=$pleiofdr_outfolder/${sm2}_vs_${sm1}_condfdr
    if [ ! -d $condfdr2 ]; then
        echo "$condfdr2 doesn't exist"
        continue
    fi

    sh $tab/fdr_fuma_snp_input.sh $conjfdr/conj.result.clump.snps.csv 0.1 0.6 $outfolder/${tag1}_vs_${tag2}_conj_snps.txt "$complex_regions"
    sh $tab/fdr_fuma_snp_input.sh $condfdr/cond.result.clump.snps.csv 0.1 0.6 $outfolder/${tag1}_vs_${tag2}_cond_snps.txt "$complex_regions"
    sh $tab/fdr_fuma_snp_input.sh $condfdr2/cond.result.clump.snps.csv 0.1 0.6 $outfolder/${tag2}_vs_${tag1}_cond_snps.txt "$complex_regions"
    cat $outfolder/${tag1}_vs_${tag2}_cond_snps.txt $outfolder/${tag2}_vs_${tag1}_cond_snps.txt | grep -v snp >> $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
    rm -f $outfolder/${tag1}_vs_${tag2}_cond_snps.txt $outfolder/${tag2}_vs_${tag1}_cond_snps.txt
    n=$((n+1))
done

if [ -s $outfolder/${tag11}_vs_${tag22}_cond_snps.txt ]; then
    echo 'snp	pval' > $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp
    cat $outfolder/${tag11}_vs_${tag22}_cond_snps.txt | cut -f1 | sort -s | uniq |  awk '{print $1,1e-10}' OFS='\t' >>  $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp
    mv $outfolder/${tag11}_vs_${tag22}_cond_snps.tmp $outfolder/${tag11}_vs_${tag22}_cond_snps.txt
fi
