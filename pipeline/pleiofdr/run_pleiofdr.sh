#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script runs pleioFDR analysis (WITHOUT ANY WARRANTY). To make it work,
# please customize values to the parameters within the "parameters to
# cumtomize" section according to your environment the first time you run it.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#

if [ $# -lt 7 ]; then
  echo "Usage:     sh run_pleiofdr.sh TRAIT1 TRAIT2 run_condfdr_flag run_conjfdr_flag run_clump_cond_flag run_clump_conj_flag run_on_cluster_flag [manh_colorlist]"
  echo "Arguments: TRAIT1 - summary statistics name of trait1"
  echo "           TRAIT2 - summary statistics name of trait2"
  echo "           run_condfdr_flag - flag to run condFDR [Y/N]"
  echo "           run_conjfdr_flag - flag to run conjFDR [Y/N]"
  echo "           run_clump_cond_flag - flag to run condFDR clumping [Y/N]"
  echo "           run_clump_conj_flag - flag to run conjFDR clumping [Y/N]"
  echo "           run_on_cluster_flag - flag to run on cluster [Y/N]"
  echo '           manh_colorlist - manhattan plot color (default "[1 0 0]")
                   [red: 1 0 0; green 0 1 0; blue: 0 0 1; orange: 1 0.5 0;
                   cyan: 0 0.75 0.75; darkgreen: 0 0.5 0; olive: 0.5 0.5 0;
                   magenta: 0.75 0 0.75]'
  echo "           manh_colorlist2 - manhattan plot color (default 1)
                   [1 3 5 7 9 11 13 15 17 19; 2 4 6 8 10 12 14 16 18 20;
                   orange sky_blue bluish_green yellow blue vermillion
                   reddish_purple black]"
  echo 'Example:   sh run_pleiofdr.sh UKB_MOOD_2019 CTG_COG_2018 N Y N Y N'
  echo 'Example:   sh run_pleiofdr.sh UKB_MOOD_2019 CTG_COG_2018 Y Y Y Y Y "[0 0 1]" 1'
  exit 0
fi

#--------------------------parameters to cumtomize------------------------#

#home directory of pleioFDR
export pleiofdr=/cluster/projects/nn9114k/yunhanc/github/pleiofdr

#home directory of python_convert
export python_convert=/cluster/projects/nn9114k/yunhanc/github/python_convert

#folder containing data of traits
export traitfolder=/cluster/projects/nn9114k/yunhanc/data/pleiofdr

#reference file for running cFDR
export ref_fdr=/cluster/projects/nn9114k/yunhanc/data/ref/ref9545380_1kgPhase3eur_LDr2p1.mat

#reference file for converting mat to csv
export ref_mat2csv=/cluster/projects/nn9114k/yunhanc/data/ref/9545380.ref

#reference for clumping results
export ref_clump=/cluster/projects/nn9114k/yunhanc/data/ref/plink_503eur

#result folder
export resultfolder=/cluster/projects/nn9114k/yunhanc/results/pleiofdr

#job account on cluster
jobaccount=nn9114k

#-------------------------------------------------------------------------#

#data of TRAIT1
export trait1file=$1.mat

#data of TRAIT2
export trait2file=$2.mat

#flag to run condFDR
export run_condfdr_flag=$3

#flag to run conjFDR
export run_conjfdr_flag=$4

#flag to run condFDR clumping 
export run_clump_cond_flag=$5

#flag to run conjFDR clumping 
export run_clump_conj_flag=$6

#flag to run on cluster
run_on_cluster_flag=$7

#manhattan plot color (matlab)
if [ $# -lt 8 ]; then
    export manh_colorlist="[1 0 0]"
else
    export manh_colorlist=$8
fi

#manhattan plot color (python)
if [ $# -lt 8 ] || [ $# -lt 9 ]; then
    export manh_colorlist2=1
else
    export manh_colorlist2=$9
fi

export legend_location="upper right"

#-------------------------------------------------------------------------#

#home directory of pipeline pleioFDR
export pipeline=$(dirname $0)

if [ "$run_on_cluster_flag" = "Y" ]; then
    sbatch --account $jobaccount $pipeline/pleiofdr.job
else
    sh $pipeline/pleiofdr.job
fi
