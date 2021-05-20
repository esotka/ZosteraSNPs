### PCA on genotype calls

rm(list=ls())
library(adegenet)
dat2.tr.num <- read.delim('data/zos.393ind.HWE.99.calls.forGenind2.txt',sep=" ")
inds <- rownames(dat2.tr.num)
pdf('output/PCA_genotypeCalls.pdf')
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

cols <- c("blue","purple","black","red")

# all
pca1 <- dudi.pca(dat2.tr.num,center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca1$li,env$site.nice,col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults + seeds",line=2)
# adults
pca2 <- dudi.pca(dat2.tr.num[env$adult.seed=="A",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca2$li,env$site.nice[env$adult.seed=="A"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("adults only",line=2)

# seeds
pca3 <- dudi.pca(dat2.tr.num[env$adult.seed=="S",],center=TRUE,scale=FALSE,scannf=FALSE)
s.class(pca3$li,env$site.nice[env$adult.seed=="S"],col=transp(cols,.6),grid=F,cstar = 0,cpoint=2)
mtext("seeds",line=1)

#library(lattice)
#print(xyplot(pca1$li[,2]~pca1$li[,1],groups=env$site.nice,auto.key=list(corner=c(1,0)),xlab="PC1",ylab="PC2"))

dev.off()


