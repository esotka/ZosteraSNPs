## analysis of spagedi SGS data
## 8 sites
rm(list=ls())
dat <- readLines("data/allSite.Spagedi.out_noClones_10categories.txt")
dat2 <- strsplit(dat,"\t")
ncat <- 10

### order of names changed to Curlew, Lynch, West, Niles
## names
dat3 <- strsplit(dat,"_")
site.depth <- paste(unlist(lapply(dat3[c(1:4,7,8,5,6)],"[[",1)),unlist(lapply(dat3[c(1:4,7,8,5,6)],"[[",2)),sep="_")

# meters
mean.dist <- c()
for (i in 3:12)
{mean.dist <- c(mean.dist,as.numeric(unlist(lapply(dat2[c(1:4,7,8,5,6)],"[[",i))))}
mean.dist <- matrix(mean.dist,nrow=ncat,ncol=8,byrow=T)

## kinship
xbar <- c()
for (i in 3:12)
{xbar <- c(xbar,as.numeric(unlist(lapply(dat2[c(1:4,7,8,5,6)+8],"[[",i))))}
xbar <- matrix(xbar,nrow=ncat,ncol=8,byrow=T)

## lower CI
CI_inf <- c()
for (i in 3:12)
{CI_inf <- c(CI_inf,as.numeric(unlist(lapply(dat2[c(1:4,7,8,5,6)+16],"[[",i))))}
CI_inf <- matrix(CI_inf,nrow=ncat,ncol=8,byrow=T)

## upper CI
CI_sup <- c()
for (i in 3:12)
{CI_sup <- c(CI_sup,as.numeric(unlist(lapply(dat2[c(1:4,7,8,5,6)+24],"[[",i))))}
CI_sup <- matrix(CI_sup,nrow=ncat,ncol=8,byrow=T)


# 2-sided test
p <- c()
for (i in 3:12)
{p <- c(p,as.numeric(unlist(lapply(dat2[c(1:4,7,8,5,6)+32],"[[",i))))}
p <- matrix(p,nrow=ncat,ncol=8,byrow=T)
p <- ifelse(p<=0.05,"*","")

pdf("output/SpagediSGS_noClones_10categories.pdf",width=10,height=4)
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1))
for(n in 1:8)
{
  plot(x=mean.dist[,n],y=xbar[,n],type="b",ylim=c(-.1,.2),xlim=c(0,45),ylab="kinship",xlab="meters",pch=20,cex=1.2)
  mtext(site.depth[n],line=-2)
  segments(-1,0,50,0,lty="dashed")
  points(x=mean.dist[,n],y=CI_inf[,n],col="red",lty="dashed",type="l")
  points(x=mean.dist[,n],y=CI_sup[,n],col="red",lty="dashed",type="l")
  text(x=mean.dist[,n],y=-0.1,p[,n],cex=1.5)
}

dev.off()


