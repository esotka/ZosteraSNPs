### AMOVA on seeds
### method 1 = pegas::vegan on mpgl numbers.
library(pegas)
rm(list=ls())
## convert beagle GL to mpgl

dat <- read.delim('data/Seeds.beagle.gz',sep="\t",header=T)
ninds <- 99
aa <- dat[,seq(1,3*ninds,3)+3]
ab <- dat[,seq(1,3*ninds,3)+4]
bb <- dat[,seq(1,3*ninds,3)+5]
tgenest <- c()
aa.prob <- aa*0
ab.prob <- ab*1
bb.prob <- bb*2
gmat <- aa.prob+ab.prob+bb.prob # row = loci; col = ind
inds <- readLines("data/ind99_seeds")
colnames(gmat) <- inds
############# ENV DATA ##############
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
env$site.nice <- factor(env$site.nice)
env$SvD <- factor(env$SvD)
env$site.nice.depth <- factor(paste(env$site.nice,env$SvD,sep="_"))
env$site.nice.depth <- factor(env$site.nice.depth,levels(env$site.nice.depth)[c(1:3,6:7,4:5)])

### only include unique genotypes (no Clones)
#noClones <- readLines("data/ind245adults_noClones_clean")
#env <- env[inds%in%noClones,]
#gmat <- gmat[,colnames(gmat)%in%noClones]

# between sites only (no deep vs shallow)

st.d <- dist(t(gmat))
#habitat <- factor(paste(env$site.nice,env$SvD))
site <- env$site.nice
print(m1 <- amova(st.d~site))
write.pegas.amova(m1,"output/amova.method1_seeds.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

#### AMOVA on genotype calls ####

library(vcfR)

dat <- read.vcfR("data/zos.393ind.HWE.99.loci19433.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
loci <- rownames(gt2)
loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
rownames(gt2) <- loci.tmp

loci.to.use <- read.table("data/loci19433")
loci.to.use$V3 <- paste(loci.to.use$V1,loci.to.use$V2,sep="_")
### subset to loci to use
gt2 <- gt2[rownames(gt2)%in%loci.to.use$V3,]
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])

colnames(geno) <- readLines("data/ind393_clean")
rownames(geno) <- rownames(gt2)

adult.seed <- substr(colnames(geno),7,7)
geno <- geno[,adult.seed=="S"]

seeds <- readLines("data/ind99_seeds")

library(parallel)

out <- mclapply(1:ncol(geno), function(i)
{
  tmp <- geno[,i]
  common.gt <- names(table(tmp))[table(tmp)==max(table(tmp))]
  tmp[is.na(tmp)] <- common.gt
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(geno)))

# Between sites

#> table(env$site.nice.depth)
#Curlew_D Curlew_S  Lynch_S   West_D   West_S  Niles_D  Niles_S 
#3       35        3        1       46        1       10 
### remove Curlew_D, Lynch_S, West_D, Niles_D

env2 <- env[env$site.nice.depth%in%c("Curlew_S","West_S","Niles_S"),]
geno2 <- geno[,colnames(geno)%in%env2$ind]

st.d <- dist(t(geno2))
#habitat <- factor(paste(env$site.nice,env$SvD))
site <- env2$site.nice
print(m1 <- amova(st.d~site))
write.pegas.amova(m1,"output/amova.method1-calls_seeds.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

## pairwise phist 
pw <- combn(c("Curlew_S","West_S","Niles_S"),2)
pw.stats.long <- list()
pw.stats.short <- c()
for (k in 1:dim(pw)[2])
{
  tmp <- geno2[,env2$site.nice.depth%in%pw[,k]]
  st.d <- dist(t(tmp))
  pops <- env2$site.nice.depth[env2$site.nice.depth%in%pw[,k]]
  m1 <- amova(st.d~pops)
  p <- m1$varcomp[[2]][1]
  pw.stats.long[[k]] <- m1; names(pw.stats.long)[k] <- paste(pw[1,k],pw[2,k],sep="-")
  sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
  phi <- getPhi(sig2)
  pw.stats.short <- rbind(pw.stats.short,data.frame(pop1=pw[1,k],pop2=pw[2,k],phiST=phi[1],p=round(p,4)))
}
print(pw.stats.long)
write.csv(pw.stats.short,"output/amova.method1-calls_seeds-pairwisePhiSt.csv")




