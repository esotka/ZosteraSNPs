#### density data - grids

dat <- read.csv("data/Permanent Plot Density Data.csv")

meta <- data.frame(
  site = c("DC","WB","NB","LP"),
  site.nice = c("Curlew","West","Niles","Lynch"))

dat$site.nice <- meta$site.nice[match(dat$Site,meta$site)]
dat$siteDepth <- paste(dat$site.nice,dat$Depth,sep="_")
library(lattice)
pdf("output/Density.pdf")
print(bwplot(Vegetative.Density~factor(siteDepth),data=dat))
print(anova(lm(Vegetative.Density~factor(Site)*factor(Depth),dat)))

print("means of vegetative height")
print(aggregate(dat$Vegetative.Density,by=list(dat$Site,dat$Depth),mean))
print(tmp$x[5:8]/tmp$x[1:4]) ### this much more

## zero truncated
print(bwplot(Flowering.Density~factor(siteDepth),data=dat))
library(lmerTest)
print(anova(lmer(Flowering.Density~Depth+(1|Site),dat)))


fl <- table(dat$siteDepth,dat$Flowering.Density>0)
fl <- fl/rowSums(fl)
p <- barplot(t(fl),xaxt="n")
mtext(c("SH","DE"),at=p,side=1)
mtext(c("Curlew","Lynch","Niles","West"),side=1,line=1.5,cex=1.5,at=c(1.3,3.7,6.1,8.5))
mtext("Presence of flowering shoot in quadrat (%)")
library(lattice)
print(histogram(~Flowering.Density | Site+Depth,data = dat,col="grey"))
dev.off()