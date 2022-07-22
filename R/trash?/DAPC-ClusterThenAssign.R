### DAPC - cluster and then assign
rm(list=ls())
library(adegenet)
#dat <- read.delim('data/Adults_noClones.beagle.gz',sep="\t",header=T)
#ninds <- 245
#aa <- dat[,seq(1,3*ninds,3)+3]
#ab <- dat[,seq(1,3*ninds,3)+4]
#bb <- dat[,seq(1,3*ninds,3)+5]
#tgenest <- c()
#aa.prob <- aa*0
#ab.prob <- ab*1
#bb.prob <- bb*2
#tot.prob <- aa.prob+ab.prob+bb.prob
#tot.prob.tr <- t(tot.prob) # no marker INFO; row = ind; col = loci

## convert beagle GL to mpgl

dat <- read.delim('data/zos.393ind.HWE.99.gl.beagle.gz',sep="\t",header=T)
ninds <- 245
aa <- dat[,seq(1,3*ninds,3)+3]
ab <- dat[,seq(1,3*ninds,3)+4]
bb <- dat[,seq(1,3*ninds,3)+5]
tgenest <- c()
aa.prob <- aa*0
ab.prob <- ab*1
bb.prob <- bb*2
gmat <- aa.prob+ab.prob+bb.prob # row = loci; col = ind
inds <- readLines("data/ind245adults_noClones_clean")
colnames(gmat) <- inds
#inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"

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
env$site.depth <- paste(env$site.nice,env$SvD,sep="_")
env$site.depth <- factor(env$site.depth)
env$site.depth <- factor(env$site.depth,levels=levels(env$site.depth)[c(1:4,7,8,5,6)])

#pdf("output/DAPC-ClusterThenAssign.pdf")
## using pops as apriori - adults
geno.a <- geno[rownames(geno)%in%env$ind,]
dapc.a <- dapc(geno.a,env$site.depth,n.pca=20,n.da=100)
cols.8 <- c("black","black","red","red","blue","blue","purple","purple")
scatter(dapc.pop,col=transp(cols.8,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
points(dapc.pop$grp.coord[,1], dapc.pop$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc.pop$grp.coord[,1], dapc.pop$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.8)

### find clusters and assign
grp <- find.clusters(tot.prob.tr) # 50 PCs and 8 clusters
cols.8diff <- c("black", #1
                        "red",#2
                        "darkgreen",#3
                        "gainsboro",#4
                        "yellow",#5
                        "deepskyblue",#6
                        "brown",#7
                        "dodgerblue4")
dapc.clust <- dapc(tot.prob.tr,grp$grp,n.pca=20,n.da=100)
scatter(dapc.clust,col=transp(cols.8diff,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
points(dapc.clust$grp.coord[,1], dapc.clust$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc.clust$grp.coord[,1], dapc.clust$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.8diff)

table.value(table(env$site.depth, grp$grp),col.lab=paste("clust",1:8))
barplot(t(table(env$site.depth,grp$grp)),col=cols.8diff)
compoplot(dapc.clust,col=cols.8diff) # 1st 50 individuals
dev.off()
