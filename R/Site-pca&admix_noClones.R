### ADMIX and PCA - no clone adults
rm(list=ls())
library(vcfR)
library(scales)
# for ADMIX
datFiles <- c("Cu.k02run1.qopt","Ly.k02run1.qopt","We.k03run1.qopt","Ni.k02run1.qopt")
indFiles <- c("Curlew66noCloneAdults","Lynch62noCloneAdults","West38noCloneAdults","Niles48noCloneAdults")
site <- c("Curlew","Lynch","West","Niles")

kcol <- c("black","lightgrey","darkgrey")
colorder <- list(
  c(1,2), #Curlew
  c(1,2), #Lynch
  c(1,2,3), #West
  c(1,2)) #Niles

# for PCA - genotype calls from bcftools


dat <- read.vcfR("data/SNPs/zos_7chr.214adults.90perc.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2

geno <- t(gt2) # row = ind; col = loci

#> dim(geno)
#[1]   214 6611
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
rownames(out2) <- rownames(geno)
colnames(out2) <- colnames(geno)

# make plots
pdf("output/Site-pca&admix_noClones.pdf",width=8,height=6)
par(mfrow=c(2,4))


for(i in 1:4)
{
  # PCA
  ind <- readLines(paste("data/ngsAdmixFiles/",indFiles[i],sep=""))
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

dat <- read.delim(paste("data/ngsAdmixFiles/",datFiles[i],sep=""),sep=" ",header = F)
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
