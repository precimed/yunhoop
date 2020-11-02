#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script extracts data from mixer analysis for supplementary bivariate
# table, but WITHOUT ANY WARRANTY.

# (c) 2020-2022 NORMENT, UiO
#-------------------------------------------------------------------------#

if [ $# -lt 1 ]; then
  echo "Usage:     sh mixer_bi_table.sh bi_csv"
  echo "Arguments: bi_csv - csv file containing bivariate info"
  exit 0
fi

bi_csv=$1
head -n1 $bi_csv | cut -f4-26,28 | sed 's/ //g'
cat $bi_csv | grep fit | cut -f4-26,28
