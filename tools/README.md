## Contents

* [identify_overlap_and_novel_loci.sh](#identify_overlap_and_novel_locish)

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
