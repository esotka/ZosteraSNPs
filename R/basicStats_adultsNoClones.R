### basic stats - adults and juveniles

## from hierfstat

rm(list=ls())
library(hierfstat)
library(vcfR)
dat = read.vcfR("data/SNPs/zos_7chr.214adults.99seeds.90perc.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])
geno <- t(geno) # row = ind; col = loci

rownames(geno) = colnames(gt2)
colnames(geno) = rownames(gt2)
#> dim(geno)
#[1]   344 6611

############# ENV DATA ##############
inds <- rownames(geno)
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

### adults no Clones

## format for hierfstat

geno.h <- data.frame(grp=env$site.depth.lifeHistory,geno)
out <- basic.stats(geno.h) # takes a minute

n <- c(table(env$site.depth.lifeHistory))
Ho <- colMeans(out$Ho,na.rm=T) # observed heterozygosities
Hs <- colMeans(out$Hs,na.rm=T) # observed gene diversities ("sometimes misleadingly called expected heterozygosity")
Fis <- colMeans(out$Fis,na.rm=T) # observed Fis ==> these were all NAs

out.summary <- data.frame(n,Ho,Hs,Fis)
out.summary$site <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",1))
out.summary$depth <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",2))
out.summary$adult.seed <- unlist(lapply(strsplit(rownames(out.summary),"_"),"[[",3))

write.csv(out.summary,"output/BasicStats_adultsNoClones.csv")

### statistical test among adults

dat <- read.csv("output/BasicStats_adultsNoClones.csv")
dat <- dat[dat$n>3,]
ad <- dat[dat$adult.seed=="A",]
anova(lm(Ho~depth,ad))
anova(lm(Hs~depth,ad))
anova(lm(Fis~depth,ad))
anova(lm(Ho~adult.seed,dat))
anova(lm(Hs~adult.seed,dat))
anova(lm(Fis~adult.seed,dat))
