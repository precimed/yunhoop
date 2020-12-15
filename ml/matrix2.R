#--------------------------- Description ----------------------------#

# Function: This script generates matrix including geno, covar, pheno data

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

# Usage:      Rscript matrix.R genodata snplist phenodata phenoname covardata outfile
#
# Arguments:  genotype  - prefix of genotype plink files
#             snplist   - a numeric or character vector indicating a subset of SNPs to be selected
#             covardata - covariate data files including one or more covariate columns
#             phenodata - phenotype data files including one or more phenotype columns
#             phenoname - a column name of phenotype in phenodata
#             outfile   - output file of analysis
#
# Example:    Rscript matrix.R testdata/test 1:1000 mostest_mental_pheno_200807.csv 20446-0.0 mostest_mental_covar_200807.csv 20446-0.0.csv

#-------------------------- Input paramters -------------------------#

# Import arguments from command line
args <- commandArgs(TRUE)

if (length(args) < 6) {
    stop("Six arguments must be supplied")
}

# genotype data
genodata = args[1]

# snp list
snplist = args[2]

# phenotype data
phenodata = args[3]

# phenotype name
phenoname = args[4]

# covariate data
covardata = args[5]

# output file
outfile = args[6]

#----------------------------- Start code ---------------------------#
source("readplink.R")
options(stringsAsFactors = FALSE)

pheno <- read.table(phenodata, header=T, sep=",", strip.white=T, as.is=T)
colnames(pheno)[1] <- 'IID'
pheno <- pheno[,c('IID',sub('^','X',sub('-','.',phenoname)))]
colnames(pheno) <- c('IID','pheno')
pheno <- pheno[!is.na(pheno$pheno),]
pheno <- pheno[pheno$pheno >= 0,]

covar <- read.table(covardata, header=T, sep=',', strip.white=T, as.is=T)
colnames(covar)[1] <- 'IID'
covar <- covar[,1:23]
covar <- na.omit(covar)

if (grepl(":", snplist, fixed = TRUE)==TRUE) {
    from = as.integer(unlist(strsplit(snplist, ":"))[1])
    to = as.integer(unlist(strsplit(snplist, ":"))[2])
    snplist = from:to
} else if (grepl(",", snplist, fixed = TRUE)==TRUE) {
    snplist = as.integer(unlist(strsplit(snplist, ",")))
} else if (grepl(";", snplist, fixed = TRUE)==TRUE) {
    snplist = as.integer(unlist(strsplit(snplist, ";")))
}

geno_snps <- read.table(paste0(genodata,'.bim'), header=F, strip.white=T, as.is=T)
geno_inds <- read.table(paste0(genodata,'.fam'), header=F, strip.white=T, as.is=T)
nsamples = nrow(geno_inds)
geno_snpStats <- get_bed_geno(genodata, snplist, nsamples)

geno <- data.frame(geno_inds[,2], geno_snpStats)
colnames(geno) <- c('IID', geno_snps$V2[snplist])
for (i in 2:ncol(geno)) {
    geno <- geno[geno[,i]!=-1,]
}

# merge dat and covar
dat <- merge(geno, covar, by="IID", all.x=F, all.y=F, sort=F)
# merge geno and pheno
dat <- merge(dat, pheno, by="IID", all.x=F, all.y=F, sort=F)
# remove rows with na
dat <- na.omit(dat)

write.table(dat, outfile, quote=F, row.names=F, col.names=T, sep=' ')
