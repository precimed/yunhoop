#!/bin/bash

#SBATCH --job-name=pleioFDR
#SBATCH --account=nn9114k
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3828M

module purge
module load MATLAB/2019a
module load Python/3.7.4-GCCcore-8.3.0
module load PLINK/1.9b_6.13-x86_64
set -o errexit # exit the scripts on any error
set -o nounset # exit the scripts on unset variables

# **************************************************************

mkdir -p $resultfolder

# reading trait1
IFS='_' read -r -a trait1_array <<< "$trait1file"
trait1name=${trait1_array[1]}

# reading trait2
IFS='_' read -r -a trait2_array <<< "$trait2file"
trait2name=${trait2_array[1]}

# configuration file
configfile=$pleiofdr/config.txt

if [ -s $pipeline/config.txt ]; then
    cp $pipeline/config.txt $configfile
else
    echo "please set correct path to the pipeline variable."
    exit 1
fi

# replace the variables
sed -i -e "/traitfolder=/c\traitfolder=$traitfolder" \
    -i -e "/manh_colorlist=/c\manh_colorlist=$manh_colorlist" \
    -i -e "/reffile=/c\reffile=$ref_fdr" $configfile

# For cond_FDR
if [ "$run_condfdr_flag" = "Y" ]; then
    # output directory
    sed -i -e "/outputdir=/c\outputdir=$resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr" \
        -i -e "/traitfile1=/c\traitfile1=$trait1file" \
        -i -e "/traitname1=/c\traitname1=$trait1name" \
        -i -e "/traitfiles=/c\traitfiles={\'$trait2file\'}" \
        -i -e "/traitnames=/c\traitnames={\'$trait2name\'}" \
        -i -e "/stattype=/c\stattype=condfdr" \
        -i -e "/fdrthresh=/c\fdrthresh=0.01" $configfile

    # run Matlab now for cond_FDR
     matlab -nodisplay -nodesktop -r "run $pleiofdr/runme.m; exit"

    # output directory
    sed -i -e "/outputdir=/c\outputdir=$resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr" \
        -i -e "/traitfile1=/c\traitfile1=$trait2file" \
        -i -e "/traitname1=/c\traitname1=$trait2name" \
        -i -e "/traitfiles=/c\traitfiles={\'$trait1file\'}" \
        -i -e "/traitnames=/c\traitnames={\'$trait1name\'}" \
        -i -e "/stattype=/c\stattype=condfdr" \
        -i -e "/fdrthresh=/c\fdrthresh=0.01" $configfile

    # run Matlab now for cond_FDR
     matlab -nodisplay -nodesktop -r "run $pleiofdr/runme.m; exit"
fi

# For conj_FDR
if [ "$run_conjfdr_flag" = "Y" ]; then
    # output directory
    sed -i -e "/outputdir=/c\outputdir=$resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr" \
        -i -e "/traitfile1=/c\traitfile1=$trait1file" \
        -i -e "/traitname1=/c\traitname1=$trait1name" \
        -i -e "/traitfiles=/c\traitfiles={\'$trait2file\'}" \
        -i -e "/traitnames=/c\traitnames={\'$trait2name\'}" \
        -i -e "/stattype=/c\stattype=conjfdr" \
        -i -e "/fdrthresh=/c\fdrthresh=0.05" $configfile

    # run Matlab now for conj_FDR
    matlab -nodisplay -nodesktop -r "run $pleiofdr/runme.m; exit"
fi

if [ "$run_clump_cond_flag" = "Y" ]; then
    # convert result.mat to result.mat.csv
    python $python_convert/fdrmat2csv.py $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/result.mat $ref_mat2csv

    # clump results
    python $python_convert/sumstats.py clump \
        --clump-field FDR \
        --force  \
        --sumstats $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/result.mat.csv \
        --bfile-chr $ref_clump/chr@ \
        --exclude-ranges 6:25119106-33854733 8:7200000-12500000 \
        --clump-p1 0.01 \
        --out $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/cond.result.clump

    # convert result.mat to result.mat.csv
    python $python_convert/fdrmat2csv.py $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/result.mat $ref_mat2csv

    # clump results
    python $python_convert/sumstats.py clump \
        --clump-field FDR \
        --force  \
        --sumstats $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/result.mat.csv \
        --bfile-chr $ref_clump/chr@ \
        --exclude-ranges 6:25119106-33854733 8:7200000-12500000 \
        --clump-p1 0.01 \
        --out $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/cond.result.clump
fi

if [ "$plot_manh_cond_flag" = "Y" ]; then
    # make manhattan plot
    python $python_convert/manhattan.py $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/result.mat.csv \
        --lead $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/cond.result.clump.lead.csv \
        --indep $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/cond.result.clump.indep.csv \
        --legend-labels "$trait1name & $trait2name" \
        --color-list $manh_colorlist2 \
        --legend-location "$legend_location" \
        --p FDR --y-label condFDR --p-thresh 0.01 \
        --out $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_condfdr/condfdr_001_manhattan_${trait1name}_${trait2name}

    # make manhattan plot
    python $python_convert/manhattan.py $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/result.mat.csv \
        --lead $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/cond.result.clump.lead.csv \
        --indep $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/cond.result.clump.indep.csv \
        --legend-labels "$trait2name & $trait1name" \
        --color-list $manh_colorlist2 \
        --legend-location "$legend_location" \
        --p FDR --y-label condFDR --p-thresh 0.01 \
        --out $resultfolder/${trait2file%%.*}_vs_${trait1file%%.*}_condfdr/condfdr_001_manhattan_${trait2name}_${trait1name}
fi

if [ "$run_clump_conj_flag" = "Y" ]; then
    # convert result.mat to result.mat.csv
    python $python_convert/fdrmat2csv.py $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/result.mat $ref_mat2csv

    # clump results
    python $python_convert/sumstats.py clump \
        --clump-field FDR \
        --force  \
        --sumstats $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/result.mat.csv \
        --bfile-chr $ref_clump/chr@ \
        --exclude-ranges 6:25119106-33854733 8:7200000-12500000 \
        --clump-p1 0.05 \
        --out $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/conj.result.clump
fi

if [ "$plot_manh_conj_flag" = "Y" ]; then
    # make manhattan plot
    python $python_convert/manhattan.py $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/result.mat.csv \
        --lead $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/conj.result.clump.lead.csv \
        --indep $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/conj.result.clump.indep.csv \
        --legend-labels "$trait1name & $trait2name" \
        --color-list $manh_colorlist2 \
        --legend-location "$legend_location" \
        --p FDR --y-label conjFDR --p-thresh 0.05 \
        --out $resultfolder/${trait1file%%.*}_vs_${trait2file%%.*}_conjfdr/conjfdr_005_manhattan_${trait1name}_${trait2name}
fi
