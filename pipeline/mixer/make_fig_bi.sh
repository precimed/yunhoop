#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script produces figures from MiXeR bivariate analysis. To make it work,
# please customize values to the parameters within the "parameters to cumtomize"
# section according to your environment the first time you run it.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh make_fig_bi.sh TRAIT1 TRAIT2"
  echo "Arguments: TRAIT1 - summary statistics name of trait1"
  echo "           TRAIT2 - summary statistics name of trait2"
  echo "Example:   sh make_fig_bi.sh UKB_MOOD_2019 CTG_COG_2018"
  exit 0
fi
#-------------------------------------------------------------------------#

## Set up job environment:
module --quiet purge   # clear any inherited modules
set -o errexit         # Exit the script on any error
set -o nounset         # Treat any unset variables as an error
module load Python/3.7.4-GCCcore-8.3.0

#--------------------------parameters to cumtomize------------------------#

mixer=/cluster/projects/nn9114k/yunhanc/github/mixer
sumstat=/cluster/projects/nn9114k/yunhanc/data/mixer
univariate_output=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/univariate
bivariate_output=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/bivariate

#-------------------------------------------------------------------------#

TRAIT1=$1
TRAIT2=$2

mkdir -p $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined

python $mixer/precimed/mixer_figures.py combine \
    --json $bivariate_output/${TRAIT1}_vs_${TRAIT2}/${TRAIT1}_vs_${TRAIT2}.fit.rep@.json \
    --out $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.fit

python $mixer/precimed/mixer_figures.py combine \
    --json $bivariate_output/${TRAIT1}_vs_${TRAIT2}/${TRAIT1}_vs_${TRAIT2}.test.rep@.json \
    --out $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.test

trait1json=$univariate_output/${TRAIT1}/combined/${TRAIT1}.test.json
trait2json=$univariate_output/${TRAIT2}/combined/${TRAIT2}.test.json

IFS='_' read -r -a trait1_array <<< ${TRAIT1}
trait1name=${trait1_array[1]}

IFS='_' read -r -a trait2_array <<< ${TRAIT2}
trait2name=${trait2_array[1]}

python $mixer/precimed/mixer_figures.py one \
    --json $trait1json $trait2json \
    --trait1 $trait1name $trait2name \
    --power-thresh 0.5 \
    --statistic mean std \
    --out $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.uni

python $mixer/precimed/mixer_figures.py two \
    --json-fit $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.fit.json \
    --json-test $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.test.json \
    --trait1 $trait1name \
    --trait2 $trait2name \
    --statistic mean std \
    --ext png svg \
    --out $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.bi

python $mixer/precimed/mixer_figures.py two \
    --json-fit $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.fit.json \
    --json-test $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.test.json \
    --trait1 $trait1name \
    --trait2 $trait2name \
    --trait1-file $sumstat/${TRAIT1}.sumstats.gz \
    --trait2-file $sumstat/${TRAIT2}.sumstats.gz \
    --statistic mean std \
    --out $bivariate_output/${TRAIT1}_vs_${TRAIT2}/combined/${TRAIT1}_vs_${TRAIT2}.bi2
