### PCA on genotype likelihoods

rm(list=ls())
library(adegenet)

## convert beagle GL to mpgl

dat <- read.delim('data/zos.393ind.HWE.99.gl.beagle.gz',sep="\t",header=T)
ninds <- 393
aa <- dat[,seq(1,3*ninds,3)+3]
ab <- dat[,seq(1,3*ninds,3)+4]
bb <- dat[,seq(1,3*ninds,3)+5]
tgenest <- c()
aa.prob <- aa*0
ab.prob <- ab*1
bb.prob <- bb*2
tot.prob <- aa.prob+ab.prob+bb.prob
tot.prob[tot.prob==0.999999] <- NA # these are uncertain assignments

### replace NAs with the average GL per locus
library(parallel)

out <- mclapply(1:ncol(tot.prob), function(i)
{
  tmp <- tot.prob[,i]
  tmp[is.na(tmp)] <- mean(tmp,na.rm=T)
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(tot.prob)))
tot.prob.tr <- t(out2) # no marker INFO; row = ind; col = loci


inds <- readLines("data/ind393_clean")
inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"

pdf('output/PCA_genotypeLikelihoods.pdf')
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
env$siteDepth <- factor(paste(env$site.nice,env$SvD,sep="_"))

cols <- c("blue","purple","black","red")

# all
pca1 <- dudi.pca(tot.prob.tr,center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca1$li,env$site.nice,col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults + seeds",line=2)
# adults
pca2 <- dudi.pca(tot.prob.tr[env$adult.seed=="A",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca2$li,env$site.nice[env$adult.seed=="A"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults only",line=2)

# adults: deep vs shallow
s.class(pca2$li,env$siteDepth[env$adult.seed=="A"],col=rep(transp(cols,.6),2),grid=F,cstar = 0,cpoint=2)

# seeds
pca3 <- dudi.pca(tot.prob.tr[env$adult.seed=="S",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca3$li,env$site.nice[env$adult.seed=="S"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("seeds",line=1)

library(lattice)
for (k in 1:4){
  env.tmp <- env[env$site.nice==levels(env$site.nice)[k] & env$adult.seed=="A",]
  env.tmp$siteDepth <- factor(env.tmp$siteDepth)
  tmp <- tot.prob.tr[inds%in%env.tmp$ind,]
  pca.tmp <- dudi.pca(tmp,center=TRUE,scale=FALSE,scannf=FALSE)
  s.class(pca.tmp$li,env.tmp$siteDepth,col=transp(c("red","black"),.6),grid=F,cstar = 0,cpoint=2)
  print(densityplot(~pca.tmp$li[,1],group=env.tmp$SvD,main=unique(env.tmp$site)))
}



dev.off()


