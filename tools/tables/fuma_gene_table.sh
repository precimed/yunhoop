#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary gene table based on FUMA
# analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 2 ]; then
  echo "Usage:     sh fuma_gene_table.sh fuma_gene_file outfile"
  echo "Arguments: fuma_gene_file - file that contains fuma gene info"
  echo "           outfile - output file"
  echo "Example:   sh fuma_gene_table.sh genes.txt trait1_vs_trait2_genes.txt"
  exit 0
fi
#-------------------------------------------------------------------------#

fuma_gene_file=$1
outfile=$2

echo "ensg	symbol	chr	start	end	minGwasP	IndSigSNPs	GenomicLocus	strand	type	entrezID	HUGO	pLI	ncRVIS	posMapSNP	posMapSNPs	posMapMaxCADD	eqtlMapSNP	eqtlMapSNPs	eqtlMapminP	eqtlMapminQ	eqtlMapts	eqtlDirection	ciMap	ciMapts" > $outfile

tail -n +2 $fuma_gene_file | awk '{if ($12>0) print $0,"Yes"; else print $0,"No"}' | awk '{if ($14>0) print $0,"Yes"; else print $0,"No"}' | awk '{print $1,$2,$3,$4,$5,$21,$22,$23,$6,$7,$8,$9,$10,$11,$24,$12,$13,$25,$14,$15,$16,$17,$18,$19,$20}' OFS='\t' >> $outfile

head -n 1 $outfile > $outfile.tmp
tail -n +2 $outfile | grep -E '(.*Yes){3}' >> $outfile.tmp
tail -n +2 $outfile | grep -v -E '(.*Yes){3}' >> $outfile.tmp

mv $outfile.tmp $outfile
