#### AMOVA on genotype calls ####
rm(list=ls())
library(pegas)
library(ade4)
library(vcfR)

dat <- read.vcfR(file="data/SNPs/zos_7chr.99seeds.90perc.recode.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
geno <- matrix(as.numeric(gt2),nrow=dim(gt2)[1])
rownames(geno) = rownames(gt2); colnames(geno) = colnames(gt2)
geno <- t(geno) # row = ind; col = loci


#> dim(geno)
#[1]  99 6611

############# ENV DATA ##############
inds2 <- rownames(geno)
env <- data.frame(ind=inds2,
                  adult.seed=substr(inds2,7,7),
                  site=substr(inds2,1,1),
                  depth=substr(inds2,2,2))
                  
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
env$site.nice.depth = factor(paste(env$site.nice,env$depth))
#> table(env$site.nice)

# Curlew  Lynch  Niles   West 
#    38      3     11     47 
#> dim(geno)
#[1]   214 3176

st.d <- dist(geno)
#habitat <- factor(paste(env$site.nice,env$SvD))
site <- env$site.nice
print(m1 <- pegas::amova(st.d~site))
write.pegas.amova(m1,"output/amova.method1-seeds.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

## pairwise phist 

pw <- combn(unique(env$site.nice),2)
pw.stats.long <- list()
pw.stats.short <- c()
for (k in 1:dim(pw)[2])
{
  tmp <- geno[env$site.nice%in%pw[,k],]
  st.d <- dist(tmp)
  pops <- env$site.nice[env$site.nice%in%pw[,k]]
  m1 <- pegas::amova(st.d~pops)
  p <- m1$varcomp[[2]][1]
  pw.stats.long[[k]] <- m1; names(pw.stats.long)[k] <- paste(pw[1,k],pw[2,k],sep="-")
  sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
  phi <- getPhi(sig2)
  pw.stats.short <- rbind(pw.stats.short,data.frame(pop1=pw[1,k],pop2=pw[2,k],phiST=phi[1],p=round(p,4)))
}
print(pw.stats.long)
write.csv(pw.stats.short,"output/amova.method1-calls_seeds-pairwisePhiSt.csv")


