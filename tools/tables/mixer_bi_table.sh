#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script tries to extract data from MiXeR analysis for supplementary
# bivariate table WITHOUT WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh mixer_bi_table.sh bi_csv"
  echo "Arguments: bi_csv - csv file containing bivariate info"
  echo "           outfile - output file"
  exit 0
fi
#-------------------------------------------------------------------------#

bi_csv=$1
outfile=$2
head -n1 $bi_csv | cut -f4-26,28 | sed 's/ //g' > $outfile
cat $bi_csv | grep fit | cut -f4-26,28 >> $outfile
