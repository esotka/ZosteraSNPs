### get genotype calls from VCF

rm(list=ls())
library(vcfR) # read.vcfR()
dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

inds <- readLines("data/ind393_clean")
loci.to.use <- read.table("data/loci19433")
loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")

### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use$V3,]







loci.per.ind <- (nrow(gt2)-colSums(is.na(gt2)))
hist(loci.per.ind)
ind.per.locus <- (343-)


gt2.noNAs <- gt2[,colSums(is.na(gt2))==0] ### get rid of all NAs

