#--------------------------- Description ----------------------------#

# Function: This script generates matrix including covar, pheno data

# Yunhan Chu (yunhanch@gmail.com)

# (c) 2020-2022 NORMENT, UiO

# Usage:      Rscript matrix.R covardata phenodata phenoname outfile
#
# Arguments:  covardata - covariate data files including one or more covariate columns
#             phenodata - phenotype data files including one or more phenotype columns
#             phenoname - a column name of phenotype in phenodata
#             outfile - output file of analysis
#
# Example:    Rscript matrix.R mostest_mental_covar_200807.csv mostest_mental_pheno_200807.csv 20446-0.0 20446-0.0.csv

#-------------------------- Input paramters -------------------------#

# Import arguments from command line
args <- commandArgs(TRUE)

if (length(args) < 4) {
    stop("Four arguments must be supplied")
}

# covariate data
covardata = args[1]

# phenotype data
phenodata = args[2]

# phenotype name
phenoname = args[3]

# output file
outfile = args[4]

#----------------------------- Start code ---------------------------#
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
#covar[,3] <- factor(covar[,3])

# merge covar and pheno
dat <- merge(covar, pheno, by="IID", all.x=F, all.y=F, sort=F)
# remove rows with na
dat <- na.omit(dat)
#rownames(dat) <- dat$Row.names
#dat <- dat[,-1]

write.table(dat, outfile, quote=F, row.names=F, col.names=T, sep=' ')
