## Contents

* [identify_novel_loci.sh](#identify_novel_locish)

## Summary of scripts

### identify_novel_loci.sh

**Function**
This script compares latest loci contained in a new file against previous
ones contained in an old file, and identifies and ouputs the novel loci.

**Usage** ``sh identify_novel_loci.sh newlocifile oldlocifile outputfile``

**Arguments**
* `newlocifile` - file that contains new loci with 4 columns (CHR, LEAD_SNP, MinBP, MaxBP)
* `oldlocifile` - file that contains old loci with 4 columns (CHR, LEAD_SNP, MinBP, MaxBP)
* `outputfile`  - output file that contains novel loci

**Example**
```
sh identify_novel_loci.sh [mood_bip_loci_cond_001.csv](https://github.com/precimed/yunhoop/blob/master/config/mood_bip_loci_cond_001.csv) [mood_gwas_loci.csv](https://github.com/precimed/yunhoop/blob/master/config/mood_gwas_loci.csv) novel_loci_cond_001_mood_bip.txt
