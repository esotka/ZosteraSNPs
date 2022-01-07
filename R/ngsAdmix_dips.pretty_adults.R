rm(list=ls())
ks <- c("k02run1","k03run1","k04run1","k05run1","k06run1","k07run1","k08run1","k09run1","k10run1")

############# ENV DATA ##############
inds <- readLines("data/ind99_seeds")
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
env$site.depth <- paste(env$site.nice,env$SvD,sep="_")
env$site.depth <- factor(env$site.depth)
env$site.depth <- factor(env$site.depth,levels=levels(env$site.depth)[c(1:3,6:7,4:5)])

kcol <- c("black", #1
          "red",#2
          "darkgreen",#3
          "gainsboro",#4
          "yellow",#5
          "deepskyblue",#6
          "brown",#7
          "dodgerblue4",#8
          "darkorchid",#9
          "burlywood")##10
colorder <- list(
  c(1,2), #ks=2
  c(1,3,2), #ks=3
  c(2,3,4,1), #ks=4
  c(1,3,5,4,2), #ks=5
  c(2,4,1,3,6,5), #ks=6
  c(4,2,3,7,6,1,5), #ks=7
  c(4,6,1,5,8,3,2,7), #ks=8
  c(6,5,3,1,7,9,2,8,4), #ks=9
  c(3,8,2,9,1,6,10,5,4,7)) #ks=10

pdf('output/ngsAdmix_dips.pretty-seeds.pdf',width=12,height=10)
par(mfrow=c(11,1),mar=c(0,0,0,0),xpd = TRUE)


for (k in 1:9)#length(ks))
{
  dat <- read.delim(paste("data/NGSadmix-seeds/",ks[k],".qopt",sep=""),sep=" ",header = F)
  dat <- dat[,-(dim(dat)[2])]
  dat <- dat[,order(colorder[[k]])]
  dat <- dat[order(env$site.depth),]
  fig <- barplot(t(dat),col=kcol[1:length(colorder[[k]])],space=0,border=NA,xlab="",ylab="",
          names.arg = rep("",nrow(dat)),horiz=F)#,main=substr(ks[k],1,3))
  mtext(substr(ks[k],1,3),side=2,line=-1.5,at=.5)
  for(i in 1:length(levels(env$site.depth)))
  {
    env.tmp <- env[order(env$site.depth),]
    x <- 1:dim(env.tmp)[1]
    xtmp <- x[env.tmp$site.depth==levels(env.tmp$site.depth)[i]]
    segments(max(xtmp),1,max(xtmp),-0.05,col="white",lwd=3)
  }
}

    fig <- barplot(t(dat),col="white",space=0,border=NA,xlab="",ylab="",
                   names.arg = rep("",nrow(dat)),horiz=F)
    for(i in 1:length(levels(env$site.depth)))
    {
      env.tmp <- env[order(env$site.depth),]
      x <- 1:dim(env.tmp)[1]
      xtmp <- x[env.tmp$site.depth==levels(env.tmp$site.depth)[i]]
      text(x=mean(xtmp),y=.5,levels(env.tmp$site.depth)[i],cex=1.5)
      segments(max(xtmp),1,max(xtmp),-0.05,col="black",lwd=3)
      
  }


dev.off()

