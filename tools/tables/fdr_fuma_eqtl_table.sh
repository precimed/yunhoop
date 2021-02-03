#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary eqtl table from FUMA analysis,
# but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh fdr_fuma_eqtl_table.sh fdr_clump_loci_file fuma_eqtl_file outfile"
  echo "Arguments: fdr_clump_loci_file - fdr clumping loci file"
  echo "           fuma_eqtl_file - file that contains fuma eqtl info"
  echo "           outfile - output file"
  echo "Example:   sh fuma_eqtl_table.sh conj.result.clump.loci.csv eqtl.txt TRAIT1_vs_TRAIT2_eqtl.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

fdr_clump_loci_file=$1
fuma_eqtl_file=$2
outfile=$3

echo 'locusnum	rsID	chr	bp	A1	A2	gene	symbol	tissue	p	signed_stats' > $outfile
cut -f1,3-8,13 $fuma_eqtl_file | tail -n +2 | sed 's/:/	/g' | awk '{if ($7==$4) print $1":"$2,$1,$2,$4,$3,$6,$11,$5,$8,$9,$10; else print $1":"$2,$1,$2,$3,$4,$6,$11,$5,$8,$9,$10}' OFS='\t' | sort -s -k1,1 > $outfile.sorted
awk '{print $2":"$4,$1,$3}' $fdr_clump_loci_file | sort -s -k1,1 > $fdr_clump_loci_file.sorted
join -1 1 -2 1 $fdr_clump_loci_file.sorted $outfile.sorted | cut -d' ' -f2-12 | sed 's/ /	/g' | sort -n -k3,3 -k4,4 >> $outfile
rm -f $fdr_clump_loci_file.sorted $outfile.sorted
