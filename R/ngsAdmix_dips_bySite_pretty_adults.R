### NGSadmix - by site

rm(list=ls())
ks <- c("Cu.k02run1.qopt","Ly.k02run1.qopt","We.k02run1.qopt","Ni.k02run1.qopt")
inds <- readLines("data/ind294adults_clean")
############# ENV DATA ##############
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice.habitat <- paste(env$site.nice,env$SvD,sep="_")
env$site.nice.habitat <- factor(env$site.nice.habitat)
levels(env$site.nice) <- levels(env$site.nice)[c(1,2,4,3)] # this is the order
pdf("output/ngsAdmix_dips_bySite_pretty-adults.pdf")
par(mfrow=c(4,1),mar=c(3,3,3,3))
for (k in 1:length(ks))
{
  
  dat <- read.delim(paste("data/NGSadmix-adults/bysite/",ks[k],sep=""),sep=" ",header = F)
  inds.tmp <- env$ind[substr(env$site.nice,1,2)==substr(ks[k],1,2)]
  hab.tmp <- substr(inds.tmp,2,2)
  site.tmp <- unique(env$site.nice[substr(env$site.nice,1,2)==substr(ks[k],1,2)])
  dat2 <- dat[order(hab.tmp),]
  dat2 <- dat2[,-(dim(dat2)[2])]
  barplot(t(dat2),col=1:nrow(dat2),space=0,border=NA,xlab="",ylab="admixture",names.arg = rep("",nrow(dat2)),main=site.tmp)
  # labels
  x <- 1:dim(dat2)[1]
  for(i in 1:2)
  {
    xtmp <- x[hab.tmp==unique(hab.tmp)[i]]
    segments(max(xtmp),1,max(xtmp),-.05,col="white",lwd=3)
    mtext(side=1,at=mean(xtmp),unique(hab.tmp)[i],cex=.75,line=1.5)
  }
}

dev.off()

