#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script runs MiXeR bivariate analysis. To make it work, please
# customize values to the parameters within the "parameters to cumtomize"
# section according to your environment the first time you run it.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh run_mixer_bi.sh TRAIT1 TRAIT2"
  echo "Arguments: TRAIT1 - summary statistics name of trait1"
  echo "           TRAIT2 - summary statistics name of trait2"
  echo "Example:   sh run_mixer_bi.sh UKB_MOOD_2019 CTG_COG_2018"
  exit 0
fi
#-------------------------------------------------------------------------#

export TRAIT1=$1
export TRAIT2=$2

#--------------------------parameters to cumtomize------------------------#

export MIXER_ROOT=/cluster/projects/nn9114k/yunhanc/github/mixer
export UNI_OUTDIR=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/univariate
export BI_OUTDIR=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/bivariate
export SUMSTAT=/cluster/projects/nn9114k/yunhanc/data/mixer
export LDFILE=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
export BIMFILE=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
export EXTRACT=/cluster/projects/nn9114k/yunhanc/data/ref/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps

#-------------------------------------------------------------------------#

sbatch $(dirname $0)/mixer_bi.job
