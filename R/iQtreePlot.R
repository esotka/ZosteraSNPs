### tree
library(ape)
library(scales)
pdf("output/iQtree.pdf")
rm(list=ls())
###############
### Adults ####
###############
tr <- read.tree('data/fasta/iqtall/iqt.contree')
plot(tr,"unrooted",show.tip.label=F)
site <- factor(substr(tr$tip.label,1,1))
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
print(levels(site)) # blue = Curlew; purple = Lynch; black = Niles; red = West 
pts <- c(19,15)[depth]
print(levels(depth)) # circle = Deep; squares = Shallow
tiplabels(cex=1,col=cols,pch=pts)
nodes.to.show <- c(266,267,268,423,424,447,269,305,306,307,308,351,352,388,400,423)#c(268,269,271,307,308,310,354,390,402,425,426,449)
nodelabels(text=tr$node.label[nodes.to.show-242],node = nodes.to.show,frame = "c",cex=.5,bg="white")
#nodelabels(frame="none")#,node = 243:(240+242))

text("Niles-Deep",x=0.12,y=0.185)
text("Niles-Shallow",x=0.02009938,y=0.185)
text("Curlew-Shallow",x=0.02849354,y=0.01)
text("Curlew-Deep",x=0.02,y=0.12)
text("Lynch-Deep",x=0.18,y=0.115)
text("Lynch-Shallow",x=0.17026156,y=0.033569727)
text("West-Deep",x=0.185,y=0.14)
text("West-Shallow",x=0.15067518,y=0.16)

###############
### Seeds ####
###############

tr <- read.tree('data/fasta/iqtall_seeds/iqt.contree')
plot(tr,"unrooted",show.tip.label=F)
site <- factor(substr(tr$tip.label,1,1))
depth <- factor(substr(tr$tip.label,2,2))
cols <- alpha(c("blue","purple","black","red"),.5)[site]
print(levels(site)) # blue = Curlew; purple = Lynch; black = Niles; red = West 
pts <- c(19,15)[depth]
print(levels(depth)) # circle = Deep; squares = Shallow
tiplabels(cex=1,col=cols,pch=pts)
nodes.to.show <- c(98,141,143)
nodelabels(text=tr$node.label[nodes.to.show-95],node = nodes.to.show,frame = "c",cex=.5,bg="white")
text("West",x=0.0,y=0.15)
text("Curlew",x=0,y=0.05)
text("Lynch",x=0.15,y=0.15)
text("Niles",x=0.17,y=0.09)
dev.off()
