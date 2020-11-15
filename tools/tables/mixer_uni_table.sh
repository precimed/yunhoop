#!/bin/bash
#--------------------------- Description ---------------------------------#
# This script extracts data from MiXeR analysis for supplementary univariate
# table, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 1 ]; then
  echo "Usage:     sh mixer_uni_table.sh uni_csv"
  echo "Arguments: uni_csv - csv file containing univariate info"
  exit 0
fi
#-------------------------------------------------------------------------#

uni_csv=$1
head -n1 $uni_csv | cut -f2-12 | sed 's/ //g'
tail -n1 $uni_csv | cut -f2-12
