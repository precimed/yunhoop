## Contents

* [mixer_uni_table.sh](#mixer_uni_tablesh)
* [mixer_bi_table.sh](#mixer_bi_tablesh)
* [fdr_fuma_snp_table.sh](#fdr_fuma_snp_tablesh)
* [fdr_fuma_loci_table.sh](#fdr_fuma_loci_tablesh)
* [fuma_genes_table.sh](#fuma_genes_tablesh)

## Summary of scripts

### mixer_uni_table.sh

**Function**
This script extracts data from mixer analysis for supplementary univariate
table.

**Usage** ``sh mixer_uni_table.sh uni_csv``

**Arguments**
* `uni_csv` - csv file containing univariate info

### mixer_bi_table.sh

**Function**
This script extracts data from mixer analysis for supplementary bivariate
table.

**Usage** ``sh mixer_bi_table.sh bi_csv``

**Arguments**
* `bi_csv` - csv file containing bivariate info

### fdr_fuma_snp_table.sh

**Function**
This script generates data for supplementary snp table (FUMA+SUMSTAT)
based on selected snp subset (fdr,r2) from cFDR/FUMA analysis.

**Usage** ``sh fdr_fuma_snp_table.sh cfdr_clump_snp_file fdr r2 fuma_snp_file sumstat1 sumstat2 outfolder``

**Arguments**
* `cfdr_clump_snp_file` - file that contains cfdr clumping snp info
* `fdr` - FDR filter for selecting snps
* `r2` - r2 filter for selecting snps
* `fuma_snp_file` - FUMA annotated snp file
* `sumstat1` - primary summary statistic file
* `sumstat2` - secondary summary statistic file
* `outfolder` - output folder

**Example**
```
sh fdr_fuma_snp_table.sh PGC_BIP_2016_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.snps.csv 0.1 0.6 snp2gene/BIP_vs_MOOD_conj_005/snps.txt sumstat/std/PGC_BIP_2016.sumstats.gz sumstat/std/UKB_MOOD_2019.sumstats.gz snp2gene
```

### fdr_fuma_loci_table.sh

**Function**
This script generates data for supplementary loci table (cFDR+FUMA+SUMSTAT)
from cFDR/FUMA analysis.

**Usage** ``sh fdr_fuma_loci_table.sh cfdr_clump_loci_file fdr_fuma_snp_table outfolder``

**Arguments**
* `cfdr_clump_loci_file` - file that contains cfdr clumping loci info
* `fdr_fuma_snp_table` - file generated by script fdr_fuma_snp_table.sh
* `outfolder` - output folder

**Example**
```
sh fdr_fuma_loci_table.sh PGC_BIP_2016_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.loci.csv snp2gene/BIP_vs_MOOD_snps_conj.txt snp2gene
```

### fuma_genes_table.sh

**Function**
This script generates data for supplementary genes table based on FUMA analysis.

**Usage** ``sh fuma_genes_table.sh fuma_genes_file outfile``

**Arguments**
* `fuma_genes_file` - file that contains fuma genes info
* `outfile` - output file

**Example**
```
sh fuma_genes_table.sh snp2gene/BIP_vs_MOOD_conj_005/genes.txt snp2gene/BIP_vs_MOOD_genes.txt
```