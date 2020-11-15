## Contents

* [run_mixer_uni.sh](#run_mixer_unish)
* [run_mixer_bi.sh](#run_mixer_bish)
* [make_fig_uni.sh](#make_fig_unish)
* [make_fig_bi.sh](#make_fig_bish)

## Summary of scripts

### run_mixer_uni.sh

**Function**
This script runs MiXeR univariate analysis. To make it work, please
customize values to the parameters within the "parameters to cumtomize"
section according to your environment the first time you run it.

**Usage** ``sh run_mixer_uni.sh TRAIT``

**Arguments**
* `TRAIT` - summary statistics name of trait

**Example**
```
sh run_mixer_uni.sh UKB_MOOD_2019
```

### run_mixer_bi.sh

**Function**
This script runs MiXeR bivariate analysis. To make it work, please
customize values to the parameters within the "parameters to cumtomize"
section according to your environment the first time you run it.

**Usage** ``sh run_mixer_bi.sh TRAIT1 TRAIT2``

**Arguments**
* `TRAIT1` - summary statistics name of trait1
* `TRAIT2` - summary statistics name of trait2

**Example**
```
sh run_mixer_bi.sh UKB_MOOD_2019 CTG_COG_2018
```

### make_fig_uni.sh

**Function**
This script produces figures from MiXeR univariate analysis. To make it work,
please customize values to the parameters within the "parameters to cumtomize"
section according to your environment the first time you run it.

**Usage** ``sh make_fig_uni.sh TRAIT``

**Arguments**
* `TRAIT` - summary statistics name of trait

**Example**
```
sh make_fig_uni.sh UKB_MOOD_2019
```

### make_fig_bi.sh

**Function**
This script produces figures from MiXeR bivariate analysis. To make it work,
please customize values to the parameters within the "parameters to cumtomize"
section according to your environment the first time you run it.

**Usage** ``sh make_fig_bi.sh TRAIT1 TRAIT2``

**Arguments**
* `TRAIT1` - summary statistics name of trait1
* `TRAIT2` - summary statistics name of trait2

**Example**
```
sh make_fig_bi.sh UKB_MOOD_2019 CTG_COG_2018
```
