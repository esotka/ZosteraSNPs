### DAPC - cluster and then assign
rm(list=ls())
library(adegenet)
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
tot.prob.tr <- t(tot.prob) # no marker INFO; row = ind; col = loci
inds <- readLines("data/ind393_clean")
inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"

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
env$site.nice <- factor(env$site.nice)

### find clusters
grp <- find.clusters(tot.prob.tr,max.n.clust=40)
dapc2 <- dapc(tot.prob.tr,grp)

### adults + seeds
cols <- c("blue","purple","black","red")
dapc2 <- dapc(tot.prob.tr,)#,grp=env$site.nice,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
mtext("adults + seeds",line=1)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols)
