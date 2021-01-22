#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary gene ontology gene-sets table
# from FUMA analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh fuma_go_table.sh fuma_gs_file outfile"
  echo "Arguments: fuma_gs_file - file that contains fuma gene-sets info"
  echo "           outfile - output file"
  echo "Example:   sh fuma_go_table.sh GS.txt trait1_vs_trait2_go.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

fuma_gs_file=$1
outfile=$2

head -n1 $fuma_gs_file > $outfile
grep '^GO_' $fuma_gs_file >> $outfile
