#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script runs MiXeR univariate analysis. To make it work, please
# customize values to the parameters within the "parameters to cumtomize"
# section according to your environment the first time you run it.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 1 ]; then
  echo "Usage:     sh run_mixer_uni.sh TRAIT"
  echo "Arguments: TRAIT - summary statistics name of trait"
  echo "Example:   sh run_mixer_uni.sh UKB_MOOD_2019"
  exit 0
fi
#-------------------------------------------------------------------------#

export TRAIT=$1

#--------------------------parameters to cumtomize------------------------#

export MIXER_ROOT=/cluster/projects/nn9114k/yunhanc/github/mixer
export OUTDIR=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/univariate
export SUMSTAT=/cluster/projects/nn9114k/yunhanc/data/mixer
export LDFILE=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
export BIMFILE=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export EXTRACT=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps

#-------------------------------------------------------------------------#

sbatch $(dirname $0)/mixer_uni.job
