## Contents

* [identify_overlap_and_novel_loci.sh](#identify_overlap_and_novel_locish)
* [check_loci_in_gwasc.sh](#check_loci_in_gwascsh)
* [identify_overlap_loci.sh](#identify_overlap_locish)

## Summary of scripts

### identify_overlap_and_novel_loci.sh

**Function**
This script tries to compare latest loci contained in a new file against previous ones contained in an old file, and identify and ouput overlapping and novel loci.

**Usage** ``sh identify_overlap_and_novel_loci.sh newlocifile oldlocifile outprefix``

**Arguments**
* `newlocifile` - file that contains new loci with 3 or 4 columns (CHR, MinBP, MaxBP, [LeadSNP])
* `oldlocifile` - file that contains old loci with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP])
* `outprefix` - prefix of the output files that contain overlapping and novel loci

**Example**
```
sh identify_overlap_and_novel_loci.sh mood_bip_loci_cond_001.csv mood_gwas_loci.csv mood_bip_cond_001

mood_bip_loci_cond_001.csv: https://github.com/precimed/yunhoop/blob/master/config/mood_bip_loci_cond_001.csv
mood_gwas_loci.csv: https://github.com/precimed/yunhoop/blob/master/config/mood_gwas_loci.csv
```

### check_loci_in_gwasc.sh

**Function**
This script tries to check whether the loci reported by fdr analysis have been registered in gwascatalog with respect to specific phenotype, in which case the loci are considered unnovel.

**Usage** ``sh check_loci_in_gwasc.sh fdr_clump_snp_file fuma_gwascatalog_file keywords``

**Arguments**
* `fdr_clump_snp_file` - fdr clumping snp file
* `fuma_gwascatalog_file` - fuma gwascatalog file
* `keywords` - keywords with respect to specific phenotype, can be multiple words combined with & and multiple patterns delimited with | within quotes
& `outfile` - output file including hits of gwascatalog

**Example**
```
sh check_loci_in_gwasc.sh conj.result.clump.snps.csv gwascatalog.txt 'high & density & lipoprotein | hdl & cholesterol' hdl_gwasc.txt'
```

### identify_overlap_loci.sh

**Function**
This script tries to identify cross-trait overlapping loci (under linux).

**Usage** ``sh identify_overlap_loci.sh list_of_loci_files outfile``

**Arguments**
* `list_of_loci_files` - file that includes paths to loci files with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP|locusnum_LeadSNP]) and tags
* `outfile` - out file that contains shared loci groups

**Example**
```
sh identify_overlap_loci.sh mood_psych_conj005_loci.txt mood_psych_conj005_shared_loci.txt

mood_psych_conj005_loci.txt: https://github.com/precimed/yunhoop/blob/master/config/mood_psych_conj005_loci.txt
mood_bip_conj005_loci.csv: https://github.com/precimed/yunhoop/blob/master/config/mood_bip_conj005_loci.csv
```
