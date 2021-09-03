rm(list=ls())
ks <- c("k02run1","k03run1","k04run1","k05run1","k06run1","k07run1","k08run1","k09run1","k10run1")

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
env$site.depth <- paste(env$site.nice,env$SvD,sep="_")

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
  c(2,1,3), #ks=3
  c(3,1,4,2), #ks=4
  c(4,1,3,2,5), #ks=5
  c(6,1,3,4,5,2), #ks=6
  c(6,1,5,2,3,4,7), #ks=7
  c(1,5,2,6,3,4,7,8), #ks=8
  c(4,8,7,1,2,3,5,9,6), #ks=9
  c(1,10,8,4,2,3,9,7,6,5)) #ks=10

pdf('output/ngsAdmix_dips.pretty-adults+seeds.pdf',width=12,height=10)
par(mfrow=c(11,1),mar=c(0,0,0,0),xpd = TRUE)


for (k in 1:9)#length(ks))
{
  dat <- read.delim(paste("data/NGSadmix-adults+seeds/",ks[k],".qopt",sep=""),sep=" ",header = F)
  dat <- dat[,-(dim(dat)[2])]
  dat <- dat[,order(colorder[[k]])]
  fig <- barplot(t(dat),col=kcol[1:length(colorder[[k]])],space=0,border=NA,xlab="",ylab="",
          names.arg = rep("",nrow(dat)),horiz=F)#,main=substr(ks[k],1,3))
  mtext(substr(ks[k],1,3),side=2,line=-1.5,at=.5)
  for(i in 1:length(unique(env$site.depth)))
  {
    x <- 1:dim(env)[1]
    xtmp <- x[env$site.depth==unique(env$site.depth)[i]]
    segments(max(xtmp),1,max(xtmp),-0.05,col="white",lwd=3)
  }
}

    fig <- barplot(t(dat),col="white",space=0,border=NA,xlab="",ylab="",
                   names.arg = rep("",nrow(dat)),horiz=F)
    for(i in 1:length(unique(env$site.depth)))
    {
      x <- 1:dim(env)[1]
      xtmp <- x[env$site.depth==unique(env$site.depth)[i]]
      text(x=mean(xtmp),y=.5,unique(env$site.depth)[i],cex=1.5)
      segments(max(xtmp),1,max(xtmp),-0.05,col="black",lwd=3)
      
  }


dev.off()

