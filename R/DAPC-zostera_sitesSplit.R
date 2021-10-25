### DAPC analysis - genotype calls - a priori groups
rm(list=ls())
library(adegenet)
pdf('output/DAPC-zostera_sitesSplit.pdf')
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
geno <- t(geno) # row = ind; col = loci

rownames(geno) <- readLines("data/ind393_clean")
colnames(geno) <- rownames(gt2)

### replace NAs with the most common genotype per locus
library(parallel)

out <- mclapply(1:ncol(geno), function(i)
{
  tmp <- geno[,i]
  common.gt <- names(table(tmp))[table(tmp)==max(table(tmp))]
  tmp[is.na(tmp)] <- common.gt
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(geno)))

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


### adults + seeds
cols <- c("blue","purple","black","red")
dapc2 <- dapc(out2,grp=env$site.nice,n.pca=20,n.da=100)
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

dapc2 <- dapc(out2,grp=env.tmp,n.pca=20,n.da=100)
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

dapc2 <- dapc(out2,grp=env.tmp,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols.tmp,0.7),pch=c(19,15),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,leg=T,lwd=3)
mtext("adults + seeds",line=1)
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols.tmp)

### adults only
adults <- out2[env$adult.seed=="A",]
adults.env <- env[env$adult.seed=="A",]
dapc3 <- dapc(adults,grp=adults.env$site.nice,n.pca=20,n.da=100)
scatter(dapc3,col=transp(cols,0.7),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=T,leg=T,lwd=3)
mtext("adults",line=1)
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc3$grp.coord[,1], dapc3$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=cols)

### adults only - coded habitats differently
adults <- out2[env$adult.seed=="A",]
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
seeds <- out2[env$adult.seed=="S",]
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
  tmp <- out2[rownames(geno)%in%env.tmp$ind,]
  grp.tmp <- factor(paste(env.tmp$SvD,env.tmp$adult.seed,sep="_"))
  dapc2 <- dapc(tmp,grp=grp.tmp,n.pca=20,n.da=100)
  scatter(dapc2,col=transp(cols[k],0.7),pch=c(19,15,21,22),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,txt.leg=c("Deep - adults","Deep - seeds", "Shallow - adults", "Shallow - seeds"),leg=T,lwd=3,main=levels(env$site.nice)[k])
  text(levels(env$site.nice)[k],x=-4,y=3)
}
  # Lynch - no deep seeds
k=2
env.tmp <- env[env$site.nice==levels(env$site.nice)[k],]
tmp <- out2[rownames(geno)%in%env.tmp$ind,]
grp.tmp <- factor(paste(env.tmp$SvD,env.tmp$adult.seed,sep="_"))
dapc2 <- dapc(tmp,grp=grp.tmp,n.pca=20,n.da=100)
scatter(dapc2,col=transp(cols[k],0.7),pch=c(19,21,22),cex=1.5,cstar=0,scree.da=F,cell=0,clab=0,solid=.4,mstree=F,txt.leg=c("Deep - adults", "Shallow - adults", "Shallow - seeds"),leg=T,lwd=3,main=levels(env$site.nice)[k])
text(levels(env$site.nice)[k],x=-4,y=3)


dev.off()
