### tree
library(ape)
library(scales)
pdf("output/iQtree.pdf",height=7,width=10)
rm(list=ls())

###############
### Seeds ####
###############

tr <- read.tree('data/fasta/iqtall_95seeds/iqt.contree')
plot(tr,"unrooted",show.tip.label=F)
site <- factor(substr(tr$tip.label,1,1))
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
print(levels(site)) # blue = Curlew; purple = Lynch; black = Niles; red = West 
pts <- c(19,17)[depth]
print(levels(depth)) # circle = Deep; triangle = Shallow
tiplabels(cex=1,col=cols,pch=pts)
#nodes.to.show <- c(102,145)
#nodelabels(text=tr$node.label[nodes.to.show-99],node = nodes.to.show,frame = "c",cex=.5,bg="white")
#text("West",x=0.02,y=0.15)
#text("Curlew",x=0,y=0.05)
#text("Lynch",x=0.15,y=0.15)
#text("Niles",x=0.15,y=0.09)


##########################
### No clones - adults ###
##########################

tr <- read.tree('data/fasta/iqtall_184adults/iqt.contree')
#tr$node.label <- as.numeric(tr$node.label)
#tr$node.label <- ifelse(tr$node.label<=100,NA,8)
plot(tr,"unrooted",show.tip.label=F)
#nodelabels(pch=tr$node.label,cex=1)
site <- factor(substr(tr$tip.label,1,1))
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
cols.unique <- alpha(c("blue","purple","black","red"),.5)
print(levels(site)) # blue = Curlew; purple = Lynch; black = 
pts <- c(19,17)[depth]
print(levels(depth)) # circle = Deep; triangle = Shallow
tiplabels(cex=1,col=cols,pch=pts)
text(x=c(0.165,0.18),y=c(0.16,0.16),c("D","S"),cex=1.5)
points(x=c(0.165,0.18),y=c(0.15,0.15),pch=c(19,17),col=cols.unique[1],cex=2)
points(x=c(0.165,0.18),y=c(0.14,0.14),pch=c(19,17),col=cols.unique[2],cex=2)
points(x=c(0.165,0.18),y=c(0.13,0.13),pch=c(19,17),col=cols.unique[3],cex=2)
points(x=c(0.165,0.18),y=c(0.12,0.12),pch=c(19,17),col=cols.unique[4],cex=2)
text(x=rep(0.19,4),y=c(0.15,0.14,0.13,0.12),c("Curlew","Lynch","Niles","West"),pos=4,col=cols.unique)
#nodelabels(frame="none",cex=.5)
#nodes.to.show <- c(216:218,334,335,353,244,318,247,281,306,282)
#nodelabels(text=tr$node.label[nodes.to.show-194],node = nodes.to.show,frame = "c",cex=.5,bg="white")
#text("Niles-Deep",x=0.16,y=0.271)
#text("Niles-Shallow",x=-0.001,y=0.260)
#text("Curlew-Shallow",x=0.012,y=0.042)
#text("Curlew-Deep",x=-0.028,y=0.174)
#text("Lynch-Deep",x=0.241,y=0.154)
#text("Lynch-Shallow",x=0.183,y=0.045)
#text("West-Deep",x=0.222,y=0.196)
#text("West-Shallow",x=0.167,y=0.211)

dev.off()
