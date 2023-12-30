#### AMOVA on genotype calls ####
rm(list=ls())
library(pegas)
library(ade4)

geno <- read.delim("data/genotypeCalls.393ind.8684loci.geno")

geno <- geno[substr(rownames(geno),7,7)=="S",]
#> dim(geno)
#[1]   245 8684

############# ENV DATA ##############
inds2 <- rownames(geno)
env <- data.frame(ind=inds2,
                  adult.seed=substr(inds2,7,7),
                  site=substr(inds2,1,1))
                  
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
#> table(env$site.nice)

#Curlew  Lynch  Niles   West 
#38      3     11     47 
# remove Lynch (n = 3)
geno = geno[!env$site.nice=="Lynch",]
#> dim(geno)
#[1]   96 8684
env = env[!env$site.nice=="Lynch",]

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


