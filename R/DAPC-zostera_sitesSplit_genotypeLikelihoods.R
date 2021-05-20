### DAPC analysis


rm(list=ls())
library(adegenet)
## convert beagle GL to mpgl

dat <- read.delim('data/zos.393ind.HWE.99.geno.beagle.gz',sep="\t",header=T)
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
inds <- readLines("data/ind393")
inds[inds=="LD3_012A.bwa.bam"] <- "LD3_12A.bwa.bam"
inds <- substr(inds,1,7)
pdf('output/DAPC-zostera_sitesSplit_genotypeLikelihoods.pdf')
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


### adults + seeds
cols <- c("blue","purple","black","red")
dapc2 <- dapc(tot.prob.tr,grp=env$site.nice,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
mtext("adults + seeds",line=1)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols)


### adults + seeds; coded depths as different symbols
env.tmp <- paste(env$site.nice,env$SvD,sep="_")
cols.tmp <- c()
for(j in 1:4){cols.tmp <- c(cols.tmp,rep(cols[j],2))}

dapc2 <- dapc(tot.prob.tr,grp=env.tmp,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols.tmp,0.7),pch=c(19,15),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
mtext("adults + seeds",line=1)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.tmp)

### adults + seeds; coded seeds as different symbols
env.tmp <- paste(env$site.nice,env$adult.seed,sep="_")
cols.tmp <- c()
for(j in 1:4){cols.tmp <- c(cols.tmp,rep(cols[j],2))}

dapc2 <- dapc(tot.prob.tr,grp=env.tmp,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols.tmp,0.7),pch=c(19,15),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
mtext("adults + seeds",line=1)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.tmp)

### adults only
adults <- tot.prob.tr[env$adult.seed=="A",]
adults.env <- env[env$adult.seed=="A",]
dapc3 <- dapc(adults,grp=adults.env$site.nice,n.pca=20,n.da=100)
scatter(dapc3,col=transp(cols,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
mtext("adults",line=1)
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols)

### adults only - coded habitats differently
adults <- tot.prob.tr[env$adult.seed=="A",]
adults.env <- env[env$adult.seed=="A",]
env.tmp <- paste(adults.env$site.nice,adults.env$SvD,sep="_")
cols.tmp <- c()
for(j in 1:4){cols.tmp <- c(cols.tmp,rep(cols[j],2))}

dapc3 <- dapc(adults,grp=env.tmp,n.pca=20,n.da=100)
scatter(dapc3,col=transp(cols.tmp,0.7),pch=c(19,15),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
mtext("adults",line=1)
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.tmp)

### seeds only
seeds <- tot.prob.tr[env$adult.seed=="S",]
seeds.env <- env[env$adult.seed=="S",]
dapc3 <- dapc(seeds,grp=seeds.env$site.nice,n.pca=20,n.da=100)
scatter(dapc3,col=transp(cols,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
mtext("seeds",line=1)
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols)


# all four combos of SvD and adult.vs.seed 
for (k in c(1,3,4)){
  env.tmp <- env[env$site.nice==levels(env$site.nice)[k],]
  tmp <- tot.prob.tr[inds%in%env.tmp$ind,]
  grp.tmp <- factor(paste(env.tmp$SvD,env.tmp$adult.seed,sep="_"))
  dapc2 <- dapc(tmp,grp=grp.tmp,n.pca=20,n.da=100)
  scatter(dapc2,col=transp(cols[k],0.7),pch=c(19,15,21,22),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,txt.leg=c("Deep - adults","Deep - seeds", "Shallow - adults", "Shallow - seeds"),leg=T,lwd=3,main=levels(env$site.nice)[k])
  text(levels(env$site.nice)[k],x=-4,y=3)
}
  # Lynch - no deep seeds
k=2
env.tmp <- env[env$site.nice==levels(env$site.nice)[k],]
tmp <- tot.prob.tr[inds%in%env.tmp$ind,]
grp.tmp <- factor(paste(env.tmp$SvD,env.tmp$adult.seed,sep="_"))
dapc2 <- dapc(tmp,grp=grp.tmp,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols[k],0.7),pch=c(19,21,22),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,txt.leg=c("Deep - adults", "Shallow - adults", "Shallow - seeds"),leg=T,lwd=3,main=levels(env$site.nice)[k])
text(levels(env$site.nice)[k],x=-4,y=3)


dev.off()
