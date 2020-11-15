## Contents

* [run_pleiofdr.sh](#run_pleiofdrsh)

## Summary of scripts

### run_pleiofdr.sh

**Function**
This script runs pleioFDR analysis. To make it work, please customize
values to the parameters within the "parameters to cumtomize" section
according to your environment the first time you run it.

**Usage** ``sh run_pleiofdr.sh TRAIT1 TRAIT2 run_condfdr_flag run_conjfdr_flag run_clump_cond_flag run_clump_conj_flag run_on_cluster_flag [manh_colorlist]``

**Arguments**
* `TRAIT1` - summary statistics name of trait1
* `TRAIT2` - summary statistics name of trait2
* `run_condfdr_flag` - flag to run condFDR [Y/N]
* `run_conjfdr_flag` - flag to run conjFDR [Y/N]
* `run_condfdr_clump_flag` - flag to run condFDR clumping [Y/N]
* `run_conjfdr_clump_flag` - flag to run conjFDR clumping [Y/N]
* `run_on_cluster_flag` - flag to run on cluster [Y/N]
* `manh_colorlist` - manhattan plot color (default red)
                   [red: 1 0 0; green 0 1 0; blue: 0 0 1; orange: 1 0.5 0;
                   cyan: 0 0.75 0.75; darkgreen: 0 0.5 0; olive: 0.5 0.5 0;
                   magenta: 0.75 0 0.75]

**Example**
```
sh run_pleiofdr.sh UKB_MOOD_2019 CTG_COG_2018 N Y N Y N
sh run_pleiofdr.sh UKB_MOOD_2019 CTG_COG_2018 Y Y Y Y Y "[0 0 1]"
```
