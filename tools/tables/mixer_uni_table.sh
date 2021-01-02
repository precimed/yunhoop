#!/bin/bash
#--------------------------- Description ---------------------------------#
# This script extracts data from MiXeR analysis for supplementary univariate
# table, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh mixer_uni_table.sh uni_csv"
  echo "Arguments: uni_csv - csv file containing univariate info"
  echo "           outfile - output file"
  exit 0
fi
#-------------------------------------------------------------------------#

uni_csv=$1
outfile=$2

head -n1 $uni_csv | cut -f2-12 | sed 's/ //g' > $outfile
tail -n +2 $uni_csv | cut -f2-12 >> $outfile
