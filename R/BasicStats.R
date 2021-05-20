## do hierfstats basic.stats


rm(list=ls())
library(adegenet)
library(hierfstat)
geno <- read.delim("data/zos.393ind.HWE.99.calls.forGenind",sep=" ",header=F)
geno <- t(geno) # row = ind; col = loci
#geno.genind <- df2genind(geno,ncode=1,NA.char=NA)
inds <- readLines("data/ind393")
inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"
inds <- substr(inds,1,7)
#names(geno.genind@all.names) <- inds
############# ENV DATA ##############
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

#### Basic stats ###
geno2 <- data.frame(paste(env$site.nice,env$adult.seed,env$SvD,sep="_"),geno)
bs <- basic.stats(geno2)


