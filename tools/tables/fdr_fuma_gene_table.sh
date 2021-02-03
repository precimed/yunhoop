#!/bin/bash
#--------------------------- Description ---------------------------------#

# This script generates data for supplementary gene table from PleioFDR/FUMA
# analysis, but WITHOUT ANY WARRANTY.

# Yunhan Chu (yunhanch@gmail.com), Guy F. L. Hindley

# (c) 2020-2022 NORMENT, UiO

#-------------------------------------------------------------------------#
if [ $# -lt 3 ]; then
  echo "Usage:     sh fdr_fuma_gene_table.sh fdr_fuma_snp_table fuma_gene_file outfile"
  echo "Arguments: fdr_fuma_snp_table - TRAIT1_vs_TRAIT2_snps.txt file generated by script fdr_fuma_snp_table.sh"
  echo "           fuma_gene_file - file that contains fuma gene info"
  echo "           outfile - output file"
  echo "Example:   sh fdr_fuma_gene_table.sh TRAIT1_vs_TRAIT2_snps.txt genes.txt TRAIT1_vs_TRAIT2_genes.txt"
  exit 0
fi
#-------------------------------------------------------------------------#
fdr_fuma_snp_table=$1
fuma_gene_file=$2
outfile=$3

rm -f $outfile
tail -n +2 $fuma_gene_file | while read line; do
    indep=`echo $line | cut -d' ' -f22`
    rm -f $outfile.tmp1 $outfile.tmp2
    echo $indep | sed 's/;/\n/g' | sed 's/:/\n/g' | sort | uniq | sort -s > $outfile.tmp1
    cut -f15 $fdr_fuma_snp_table | tail -n +2 | sort | uniq | sort -s > $outfile.tmp2
    if [ `join $outfile.tmp1 $outfile.tmp2 | wc -l` -gt 0 ]; then
        echo $line >> $outfile
    fi
done

cat $outfile | awk '{if($12>0) print $0,"Yes"; else print $0,"No"}' | awk '{if($14>0) print $0,"Yes"; else print $0,"No"}' | awk '{print $1,$2,$3,$4,$5,$21,$22,$23,$6,$7,$8,$9,$10,$11,$24,$12,$13,$25,$14,$15,$16,$17,$18,$19,$20}' | awk '{if($15=="Yes" && $18=="Yes" && $24=="Yes") print $0,"Yes"; else print $0,"No"}' > $outfile.tmp

echo "ensg	symbol	chr	start	end	minGwasP	IndSigSNPs	GenomicLocus	strand	type	entrezID	HUGO	pLI	ncRVIS	posMapSNP	posMapSNPs	posMapMaxCADD	eqtlMapSNP	eqtlMapSNPs	eqtlMapminP	eqtlMapminQ	eqtlMapts	eqtlDirection	ciMap	ciMapts	credible" > $outfile

awk '$26=="Yes"' $outfile.tmp | sed 's/ /	/g' >> $outfile
awk '$26=="No"' $outfile.tmp | sed 's/ /	/g' >> $outfile
rm -f $outfile.tmp $outfile.tmp1 $outfile.tmp2
