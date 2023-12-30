### tree
library(ape)
library(scales)
pdf("output/iQtree.pdf",height=7,width=10)
rm(list=ls())


par(mfrow=c(1,2),mar=c(0,0,0,0))

##########################
### No clones - adults ###
##########################

tr <- read.tree('data/fasta/iqtall_184adults/iqt.contree')
plot(tr,"unrooted",show.tip.label=F)
site <- factor(substr(tr$tip.label,1,1)); site = factor(site,levels=levels(site)[c(1,2,4,3)])
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
cols.unique <- alpha(c("blue","purple","black","red"),.5)
print(levels(site)) # blue = Curlew; purple = Lynch; black = 
pts <- c(19,17)[depth]
print(levels(depth)) # circle = Deep; triangle = Shallow
tiplabels(cex=2,col=cols,pch=pts)
text(x=0.07,y=0.16,"Adults",cex=3,pos=3)

###############
### Seeds ####
###############
tr <- read.tree('data/fasta/iqtall_95seeds/iqt.contree')
plot(tr,"unrooted",show.tip.label=F)
site <- factor(substr(tr$tip.label,1,1)); site = factor(site,levels=levels(site)[c(1,2,4,3)])
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
print(levels(site)) # blue = Curlew; purple = Lynch; black = Niles; red = West 
pts <- c(19,17)[depth]
print(levels(depth)) # circle = Deep; triangle = Shallow
tiplabels(cex=2,col=cols,pch=pts)
text(x=0.12,y=0.21,"Seeds",cex=3,pos=3)

## legend ##
text(x=c(0.11,0.13),y=c(0.08,0.08),c("D","S"),cex=1.5)
points(x=c(0.11,0.13),y=c(0.065,0.065),pch=c(19,17),col=cols.unique[1],cex=2)
points(x=c(0.11,0.13),y=c(0.05,0.05),pch=c(19,17),col=cols.unique[2],cex=2)
points(x=c(0.11,0.13),y=c(0.035,0.035),pch=c(19,17),col=cols.unique[3],cex=2)
points(x=c(0.11,0.13),y=c(0.02,0.02),pch=c(19,17),col=cols.unique[4],cex=2)
text(x=rep(0.14,4),y=c(0.065,0.05,0.035,0.02),c("Curlew","Lynch","West","Niles"),pos=4,col=cols.unique)
rect(0.10,0.015,0.165,0.085)
dev.off()
