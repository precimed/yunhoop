#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script produces figures from MiXeR univariate analysis. To make it work,
# please customize values to the parameters within the "parameters to cumtomize"
# section according to your environment the first time you run it.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 1 ]; then
  echo "Usage:     sh make_fig_uni.sh TRAIT"
  echo "Arguments: TRAIT - summary statistics name of trait"
  echo "Example:   sh make_fig_uni.sh UKB_MOOD_2019"
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
univariate_output=/cluster/projects/nn9114k/yunhanc/results/mixer_1.3/univariate

#-------------------------------------------------------------------------#

TRAIT=$1

mkdir -p $univariate_output/$TRAIT/combined

python $mixer/precimed/mixer_figures.py combine \
    --json $univariate_output/$TRAIT/$TRAIT.fit.rep@.json \
    --out $univariate_output/$TRAIT/combined/$TRAIT.fit

python $mixer/precimed/mixer_figures.py combine \
    --json $univariate_output/$TRAIT/$TRAIT.test.rep@.json \
    --out $univariate_output/$TRAIT/combined/$TRAIT.test

IFS='_' read -r -a trait_array <<< $TRAIT
traitname=${trait_array[1]}

python $mixer/precimed/mixer_figures.py one \
    --json $univariate_output/$TRAIT/combined/$TRAIT.fit.json \
    --trait1 $traitname \
    --statistic mean std \
    --out $univariate_output/$TRAIT/combined/$TRAIT.fit

python $mixer/precimed/mixer_figures.py one \
    --json $univariate_output/$TRAIT/combined/$TRAIT.test.json \
    --trait1 $traitname \
    --statistic mean std \
    --power-thresh 0.5 \
    --out $univariate_output/$TRAIT/combined/$TRAIT
