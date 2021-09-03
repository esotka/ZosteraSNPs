### Are there clones?

rm(list=ls())
library(adegenet)
library(poppr)
geno <- read.table('data/zos.393ind.HWE.99.calls.forGenind',sep=" ")
geno <- t(geno)
colnames(geno) <- readLines("data/loci19433")
rownames(geno) <- readLines("data/ind393_clean")

geno.genind <- df2genind(geno,ncode=1,NA.char=NA)
names(geno.genind@all.names) <- colnames(geno)

### adults only ###
inds <- readLines("data/ind393_clean")
geno.genind.adult <- geno.genind[substr(inds,7,7)=="A"]

### meta of adults####
inds.adult <- inds[substr(inds,7,7)=="A"]
env <- data.frame(ind=substr(inds.adult,1,7),
                  adult.seed=substr(inds.adult,7,7),
                  site=substr(inds.adult,1,1),
                  SvD=tolower(substr(inds.adult,2,2)),
                  grid=substr(inds.adult,3,3),
                  ind2=substr(inds.adult,5,6))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "cur"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "nil"
env$site.nice[env$site.nice=="W"] <- "wes"
env$site.nice[env$site.nice=="L"] <- "lyn"
env$site.nice <- factor(env$site.nice)
env$quad.name <- paste(env$site.nice,".",env$SvD,env$grid,"_",env$ind2,sep="")
env$site.nice.depth <- paste(env$site.nice,env$SvD,"_")

### number of genotypes across populations
pop(geno.genind.adult) <- env$site.nice.depth
mlg(geno.genind.adult)

#############################
# Number of Individuals:  294 
# Number of MLG:  294 
#############################



