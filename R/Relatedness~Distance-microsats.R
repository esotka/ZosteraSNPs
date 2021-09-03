### microsatellite genetic distances ~ geographic distances

### enter .csv file into genind.
library(poppr)
rm(list=ls())
dat <- read.csv("data/DC_grids_final0521_microsats.csv")
rownames(dat) <- dat$sample; dat <- dat[,-1]
locusname <- colnames(dat)[seq(1,19,2)]
allele1 <- seq(1,19,2); allele2 <- seq(2,20,2)

### metadata genetics
meta <- data.frame(ind=rownames(dat),
                   site=substr(rownames(dat),1,2),
                   site.nice="cur",
                   depth=tolower(substr(rownames(dat),4,4)),
                   grid=substr(rownames(dat),5,5),
                   quad=substr(rownames(dat),7,8))
meta$quad <- ifelse(nchar(meta$quad)==1,paste("0",meta$quad,sep=""),meta$quad)
meta$quad.name <- paste(meta$site.nice,".",meta$depth,meta$grid,"_",meta$quad,sep="")


out <- c()
for (i in 1:length(locusname))
{
tmp1 <- as.character(dat[,allele1[i]])
tmp1 <- ifelse(nchar(tmp1)==3,tmp1,paste("0",tmp1,sep=""))
tmp2 <- as.character(dat[,allele2[i]])
tmp2 <- ifelse(nchar(tmp2)==3,tmp2,paste("0",tmp2,sep=""))
out <- cbind(out,paste(tmp1,tmp2,sep=""))
}
rownames(out) <- rownames(dat)
usat <- df2genind(out,ncode=3)



#### get geographic distances for these sites ####

m <- read.csv("output/GeographicDistanceIndiv.csv")
rownames(m) <- m[,1]
m <- m[,-1]
m2 <- m[rownames(m)%in%meta$quad.name,colnames(m)%in%meta$quad.name]
m3 <- m2[lower.tri(m2)]

#### Put geno into the same order as geo dist ###
#### calculate relatedness ###

usat2 <- usat[match(rownames(m2),meta$quad.name)]
nei <- nei.dist(usat)
edw <- edwards.dist(usat)
rey <- reynolds.dist(usat)
### make labels
tmp <- matrix(rownames(m2),89,89)
tmp1 <- tmp[lower.tri(tmp)]
tmp2 <- t(tmp)
tmp2 <- tmp2[lower.tri(tmp2)]

#### analysis ###
out <- data.frame(nei=as.vector(nei),
                  edw=as.vector(edw),
                  rey=as.vector(rey),
                  m=m3,ind1=tmp1,ind2=tmp2)
out$depth <- ifelse(substr(out$ind1,5,5)==substr(out$ind2,5,5),substr(out$ind1,5,5),"NA") 
out <- out[!out$depth=="NA",]


library(lattice)
pdf("output/Relatedness~Distance-microsats.pdf")
dist.names <- c("nei","edw","rey")
for(i in 1:3)
{
  print(xyplot(out[,dist.names[i]]~m,type=c("p","r"),data=out,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
  print(xyplot(out[,dist.names[i]]~m | depth,type=c("p","r"),data=out,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
}
dev.off()

### stats ####
print("All samples")
out.stats <- c()
for (j in 1:2)
{
  for(i in 1:3)
  {
    tmp <- data.frame(x=out$m[out$depth==unique(out$depth)[j]]
                      ,y=out[out$depth==unique(out$depth)[j],dist.names[i]])
    tmp.stats <- cor.test(tmp$x,tmp$y,method="kendall")
    out.stats <- rbind(out.stats,
                       data.frame(depth=unique(out$depth)[j],
                                  dist.names=dist.names[i],
                                  n=dim(tmp)[1],
                                  stat=tmp.stats$statistic,
                                  p=round(tmp.stats$p.value,3),
                                  signif=ifelse(tmp.stats$p.value<0.05,"*","")))
  }
}
print(out.stats)

### pull out clones ###

out2 <- out[out$nei>0.01,]
library(lattice)
pdf("output/Relatedness~Distance-microsats-noClones.pdf")
dist.names <- c("nei","edw","rey")
for(i in 1:3)
{
  print(xyplot(out2[,dist.names[i]]~m,type=c("p","r"),data=out2,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
  print(xyplot(out2[,dist.names[i]]~m | depth,type=c("p","r"),data=out2,ylab="genetic distance",main=dist.names[i],cex=0.5,pch=20,lwd=4,col="grey"))
}
dev.off()

### stats ####
print("Removed clones")
out2.stats <- c()
for (j in 1:2)
{
  for(i in 1:3)
  {
    tmp <- data.frame(x=out2$m[out2$depth==unique(out2$depth)[j]]
                      ,y=out2[out2$depth==unique(out2$depth)[j],dist.names[i]])
    tmp.stats <- cor.test(tmp$x,tmp$y,method="kendall")
    out2.stats <- rbind(out2.stats,
                       data.frame(depth=unique(out2$depth)[j],
                                  dist.names=dist.names[i],
                                  n=dim(tmp)[1],
                                  stat=tmp.stats$statistic,
                                  p=round(tmp.stats$p.value,3),
                                  signif=ifelse(tmp.stats$p.value<0.05,"*","")))
  }
}
print(out2.stats)


