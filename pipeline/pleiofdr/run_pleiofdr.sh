jobaccount=nn9114k

#home directory of pleiofdr
export pleiofdr=/cluster/projects/nn9114k/yunhanc/github/pleiofdr

#home directory of python_convert
export python_convert=/cluster/projects/nn9114k/yunhanc/github/python_convert

#folder containing data of traits
export traitfolder=/cluster/projects/nn9114k/yunhanc/data/pleiofdr

#data of trait1
export trait1file=CTG_COG_2018.mat

#data of trait2
export trait2file=UKB_MOOD_2019.mat

#flag to run condfdr
export run_condfdr_flag='Y'

#flag to run conjfdr
export run_conjfdr_flag='Y'

#flag to run clumping 
export run_clump_flag='Y'

#reference file for running fdr
export ref_fdr=/cluster/projects/nn9114k/yunhanc/data/ref/ref9545380_1kgPhase3eur_LDr2p1.mat

#reference file for converting mat to csv
export ref_mat2csv=/cluster/projects/nn9114k/yunhanc/data/ref/9545380.ref

#reference for clumping results
export ref_clump=/cluster/projects/nn9114k/yunhanc/data/ref/plink_503eur

#result folder
export resultfolder=/cluster/projects/nn9114k/yunhanc/results/pleiofdr

sh pleiofdr.job
#sbatch --account $jobaccount pleiofdr.job
