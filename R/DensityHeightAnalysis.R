#####################################################
## Analysis of shoot morphology in permanent quadrats 
## 25 cm x 25 cm quadrats set at every other node 
## spatially paired with our genetic sampling 
#####################################################

rm(list=ls())
dat.all <- read.csv("data/perm.quad.density.2019.csv")
# parse to include just the June data from the grids; the August data was collected when we collected flowers for the outplant experiment later that year; not spatially paired with SNPs data, but same rough location 

dat <- dat.all[dat.all$month=="June",]
dat$biomass <- as.numeric(dat$biomass)
meta <- data.frame(
  site = c("DC","WB","NB","LP"),
  site.nice = c("Curlew","West","Niles","Lynch"))
dat$site <- factor(dat$site)

dat$site.nice <- meta$site.nice[match(dat$site,meta$site)]
dat$siteDepth <- paste(dat$site.nice,dat$depth,sep="_")
library(lattice)
pdf("output/Density.pdf")
print(histogram(~veg.density | site+depth,data = dat,col="grey",breaks=20))
print(bwplot(veg.density~factor(siteDepth),data=dat))
print(anova(lm(veg.density~factor(site)*factor(depth),dat)))

print("means of vegetative density")
print(aggregate(dat$veg.density,by=list(dat$site.nice,dat$depth),mean))

## zero truncated
print(bwplot(Flowering.Density~factor(siteDepth),data=dat))
library(lmerTest)
print(anova(lmer(Flowering.Density~Depth*Site+(1|Permanent.Quadrat),dat)))


fl <- table(dat$siteDepth,dat$Flowering.Density>0)
fl <- fl/rowSums(fl)
p <- barplot(t(fl),xaxt="n")
mtext(c("SH","DE"),at=p,side=1)
mtext(c("Curlew","Lynch","Niles","West"),side=1,line=1.5,cex=1.5,at=c(1.3,3.7,6.1,8.5))
mtext("Presence of flowering shoot in quadrat (%)")
print(histogram(~log(Flowering.Density) | Site+Depth,data = dat,col="grey"))
dev.off()
