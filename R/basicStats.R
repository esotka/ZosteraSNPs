### basic stats - adults and juveniles

## from hierfstat

rm(list=ls())
library(hierfstat)
geno <- read.table('data/zos.393ind.HWE.99.calls.forGenind',sep=" ")
geno <- t(geno) # row = ind; col = loci
colnames(geno) <- readLines("data/loci19433")


############# ENV DATA ##############
inds <- readLines("data/ind393_clean")
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.depth.lifeHistory <- paste(env$site.nice,env$SvD,env$adult.seed,sep="_")

## format for hierfstat

geno.h <- data.frame(grp=env$site.depth.lifeHistory,geno)
out <- basic.stats(geno.h) # takes a minute

n <- colMeans(out$n.ind.samp)
Ho <- colMeans(out$Ho) # observed heterozygosities
Hs <- colMeans(out$Hs) # observed gene diversities ("sometimes misleadingly called expected heterozygosity")
Fis <- colMeans(out$Fis) # observed Fis ==> these were all NAs

out.summary <- data.frame(n,Ho,Hs,Fis)
out.summary$site <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",1))
out.summary$depth <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",2))
out.summary$adult.seed <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",3))

write.csv(out.summary,"output/BasicStats.csv")

### statistical test among adults

dat <- read.csv("output/BasicStats.csv")
dat <- dat[dat$n>3,]
ad <- dat[dat$adult.seed=="A",]
anova(lm(Ho~depth,ad))
anova(lm(Hs~depth,ad))
anova(lm(Ho~adult.seed,dat))
