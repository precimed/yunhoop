module purge
module load Python/3.8.2-GCCcore-9.3.0

phenoname=20446-0.0
pheno=/cluster/projects/p33/projects/mental/pheno/mostest_mental_pheno_200807.csv
covar=/cluster/projects/p33/projects/mental/covar/mostest_mental_covar_200807.csv
geno=/cluster/projects/p33/projects/mental/geno_hapmap/all_in_one/ukb_imp_v3_qc
gwas=/tsd/p33/data/durable/groups/biostat/mental/grant-erc/hapmap-raw/20446-0.0.logistic.csv
outfolder=/tsd/p33/data/durable/characters/yunhanc/mental

snplist=`cut -f1,4 $gwas | awk '$2<1e-6 {print NR}' | sed ':a;N;$!ba;s/\n/,/g'`
echo $snplist

#Rscript matrix.R $covar $pheno $phenoname $outfolder/$phenoname.csv
#python3 model.py $outfolder/$phenoname.csv
Rscript matrix2.R $geno $snplist $pheno $phenoname $covar $outfolder/$phenoname-2.csv
python3 model.py $outfolder/$phenoname-2.csv
