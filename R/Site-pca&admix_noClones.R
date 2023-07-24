### ADMIX and PCA - no clone adults
rm(list=ls())
library(vcfR)
library(scales)
# for ADMIX
datFiles <- c("Cu.k03run1.qopt","Ly.k03run1.qopt","We.k02run1.qopt","Ni.k03run1.qopt")
indFiles <- c("CurlewAdults.83inds","LynchAdults.87inds","WestAdults.55inds","NilesAdults.69inds")
site <- c("Curlew","Lynch","West","Niles")

kcol <- c("black","lightgrey","darkgrey")
colorder <- list(
  c(1,2,3), #Curlew
  c(1,2,3), #Lynch
  c(2,1), #West
  c(1,2,3)) #Niles

# for PCA - genotype calls from bcftools


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
#> dim(geno)
#[1]   393 19432
### remove clones ###

noClones <- readLines("data/ind245adults_noClones_clean")
geno <- geno[rownames(geno)%in%noClones,]
#> dim(geno)
#[1]   245 19432

#### only include loci that had NAs in <50% of individuals
prop.ind <- colSums(is.na(geno))/nrow(geno)
geno <- geno[,prop.ind<=0.5]
#> dim(geno)
#[1]   245 12704

### replace NAs with the most common genotype per locus
library(parallel)

out <- mclapply(1:ncol(geno), function(i)
{
  tmp <- geno[,i]
  common.gt <- names(table(tmp))[table(tmp)==max(table(tmp))]
  tmp[is.na(tmp)] <- common.gt
  tmp
},mc.cores = 4)

out2 <- as.data.frame(matrix(as.numeric(unlist(out)),nrow=nrow(geno))) # rows = 393 ind cols = loci
rownames(out2) <- noClones
colnames(out2) <- colnames(geno)

# make plots
pdf("output/Site-pca&admix_noClones.pdf",width=8,height=6)
par(mfrow=c(2,4))


for(i in 1:4)
{
  # PCA
  ind <- readLines(paste("data/ngsadmix-noClones-site/",indFiles[i],sep=""))
  tmp <- out2[rownames(out2)%in%ind,] # rows = pop; cols=loci
  depth <- substr(rownames(tmp),2,2); depth <- factor(depth)
  pc.cr <- prcomp(tmp)
  col.sub <- c("black","red")[depth]
  #print(depth)
  par(mar=c(3,3,4,2))
  plot(pc.cr$x[,1],pc.cr$x[,2],cex=0,xaxt="none",yaxt="none",ylab="",xlab="")
  mtext(side=1,line=1,"PC1"); mtext(side=2,line=1,"PC2")
  points(pc.cr$x[,1],pc.cr$x[,2],pch=20,cex=1.5,col=alpha(col.sub,.5))
  for (j in 1:2){
    xbar.x <- mean(pc.cr$x[depth==levels(depth)[j],1])
    xbar.y <- mean(pc.cr$x[depth==levels(depth)[j],2])
    tmp.pc.cr <- pc.cr$x[depth==levels(depth)[j],]
    segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp.pc.cr[,1],y1 = tmp.pc.cr[,2],
             col=c("black","red")[j],lwd=.5)
  }
  
  # ADMIX

dat <- read.delim(paste("data/ngsadmix-noClones-site/",datFiles[i],sep=""),sep=" ",header = F)
depth <- substr(ind,2,2); depth <- factor(depth)
dat <- dat[,-(dim(dat)[2])]
dat <- dat[,order(colorder[[i]])]
par(mar=c(2,2,2,2))
fig <- barplot(t(dat),col=kcol[1:length(colorder[[i]])],space=0,border=NA,xlab="",ylab="Proportion",cex.lab=1.2,names.arg = rep("",nrow(dat)),horiz=F,ylim=c(-.1,1.1))
## annotations
titles = paste(c("A.","B.","C.","D."),site)
mtext(titles[i],at=-15,cex=1.5)
x <- 1:dim(dat)[1]
x1 <- x[depth==levels(depth)[1]]
x2 <- x[depth==levels(depth)[2]]
segments(c(0,max(x1),max(x2)),0,c(0,max(x1),max(x2)),-.3,lwd=1,col="black")
mtext(side=1,at=mean(x1),levels(depth)[1],cex=1,line=-1)
mtext(side=1,at=mean(x2),levels(depth)[2],cex=1,line=-1)

}

dev.off()
