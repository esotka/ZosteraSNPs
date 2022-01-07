##### microsat clones on map
library(scales) # brewer
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
meta$site.nice.depth <- paste(meta$site.nice,meta$depth,sep="_")
meta$site.nice.depth <- factor(meta$site.nice.depth)


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
usat@pop <- meta$site.nice.depth

#usat2 <- usat[match(rownames(m2),meta$quad.name)]

# WHERE ARE THE USAT CLONES?
# here are geographic positions
latlon <- read.csv("output/IndivLatLon.csv")
tmp <- ifelse(latlon$ind<10,paste("0",latlon$ind,sep=""),latlon$ind)
latlon$ind <- tmp
latlon$quad.name <- paste(latlon$grid,latlon$ind,sep = "_")

mlgID <- list()
threshold = 0.03
pdf("output/usat_ClonesOnMap=threshold=0.03.pdf",width=10,height=4)
#quartz(width=3,height=4)
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1), c(2,1,0,0))
for (i in 1:2)
{
  usat.sub <- usat[usat@pop==levels(meta$site.nice.depth)[i]]
  bwd <- bitwise.dist(usat.sub,missing_match=T)
  hist(bwd,col="grey",breaks=50,main=levels(meta$site.nice.depth)[i])
  segments(threshold,-1,threshold,100,col="red")
  # Visual inspection suggests threshold of 0.01
  mlgID[[i]] <- mlg.filter(usat.sub,threshold = threshold,distance=bwd)
}

for (i in 1:2)
{
  cloneID <- mlgID[[i]]
  tmpID <- meta$ind[meta$site.nice.depth==levels(meta$site.nice.depth)[i]]
  tmpQuadName <- meta$quad.name[meta$site.nice.depth==levels(meta$site.nice.depth)[i]]
  tmpLat <- latlon$lat[match(tmpQuadName,latlon$quad.name)]
  tmpLon <- latlon$lon[match(tmpQuadName,latlon$quad.name)]
  all <- data.frame(cloneID,tmpID,tmpQuadName,tmpLat,tmpLon)
  plot(tmpLon,tmpLat,main=levels(meta$site.nice.depth)[i])
  ### colored are clones
  cloneID.plot <- names(table(cloneID))[table(cloneID)>1]
  cols <- brewer.pal(9,"YlGnBu")
  print(cloneID.plot)
    for (k in 1:length(cloneID.plot))
    {
      all.sub <- all[all$cloneID==cloneID.plot[k],]
      x <- combn(all.sub$tmpLon,2)
      y <- combn(all.sub$tmpLat,2)
      if (dim(all.sub)[1]==2) {
        segments(all.sub$tmpLon[1],all.sub$tmpLat[1],all.sub$tmpLon[2],all.sub$tmpLat[2],col=cols[k],lwd=6)}
      else {
        for (kk in 1:dim(all.sub)[1])
        {
      segments(x[1,kk],y[1,kk],x[2,kk],y[2,kk],col=cols[k],lwd=6)}
      }
    }
}

dev.off()  
