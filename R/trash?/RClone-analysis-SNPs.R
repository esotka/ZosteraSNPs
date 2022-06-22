### Using RClone to estimate MLG and MLL
rm(list=ls())

library(RClone)
library(adegenet)
library(poppr)
library(adegenet)
library(vcfR)

dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

loci.to.use <- read.table("data/loci19433")
loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")
### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use$V3,]
geno <- matrix(gt2,nrow=dim(gt2)[1])
geno <- t(geno) # row = ind; col = loci

## adults only
inds <- readLines("data/ind393_clean")
env <- data.frame(ind=inds,
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
env$SvD <- factor(env$SvD)

adults <- geno[env$adult.seed=="A",]
for(i in 1:dim(adults)[2])
{
  adults[,i] <- as.character(adults[,i])
}
# convert to RClone format

forRClone <- convert_GC(adults[,1:100],1,"/")
