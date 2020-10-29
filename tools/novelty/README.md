## Contents

* [identify_overlap_and_novel_loci.sh](#identify_overlap_and_novel_locish)
* [check_loci_in_gwasc.sh](#check_loci_in_gwascsh)

## Summary of scripts

### identify_overlap_and_novel_loci.sh

**Function**
This script compares latest loci contained in a new file against previous
ones contained in an old file, and identifies and ouputs the overlapping
and novel loci.

**Usage** ``sh identify_overlap_and_novel_loci.sh newlocifile oldlocifile outprefix``

**Arguments**
* `newlocifile` - file that contains new loci with 3 or 4 columns (CHR, MinBP, MaxBP, [LeadSNP])
* `oldlocifile` - file that contains old loci with 3 or 4 columns (CHR, MinBP, MaxBP, [leadSNP])
* `outprefix`   - prefix of the output files that contain overlapping and novel loci

**Example**
```
sh identify_overlap_and_novel_loci.sh mood_bip_loci_cond_001.csv mood_gwas_loci.csv mood_bip_cond_001

mood_bip_loci_cond_001.csv: https://github.com/precimed/yunhoop/blob/master/config/mood_bip_loci_cond_001.csv
mood_gwas_loci.csv: https://github.com/precimed/yunhoop/blob/master/config/mood_gwas_loci.csv
```

### check_loci_in_gwasc.sh

**Function**
This script checks whether the loci reported by cfdr analysis have been
registered in gwascatalog with respect to specific phenotype, if so the
locus is considered unnovel.

**Usage** ``sh check_loci_in_gwasc.sh cfdr_clump_snps_file fuma_gwascatalog_file keyword``

**Arguments**
* `cfdr_clump_snps_file` - file that contains snp info from cfdr clumping
* `fuma_gwascatalog_file` - file that contains gwascatalog info from FUMA
* `keyword` - keyword such as name of part of the name of investigated trait

**Example**
```
sh check_loci_in_gwasc.sh PGC_BIP_2016_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.snps.csv snp2gene/BIP_vs_MOOD_conj_005/gwascatalog.txt bipolar

sh check_loci_in_gwasc.sh PGC_MDD_2018_with23andMe_noUKBB_vs_UKB_MOOD_2019_conjfdr/conj.result.clump.snps.csv snp2gene/MDD_vs_MOOD_conj_005/gwascatalog.txt depress
```
