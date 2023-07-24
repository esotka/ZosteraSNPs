### field data vs Nb (from SPAGEDi) and Ne (from NeEstimator of adults)
rm(list=ls())
dat = read.table("data/Nb_FieldDensity.txt",header=T)
library(lattice)
print(xyplot(Ne_SPAGEDI~FlowerPlantDensityQuad | SvD,dat,type=c("p","r")))
print(xyplot(Ne_SPAGEDI~PlantDensityQuad | SvD,dat,type=c("p","r")))
print(xyplot(Ne_SPAGEDI~clonalDiversity ,dat,type=c("p","r")))

size = 40*40 #= 1600 m2 ### spatial spread 
quad = 0.0625 # per quadrat m2
den = size/quad # = 1600/0.0625 # num. of quadrats per 40 x 40m space 
dat$Nc = dat$PlantDensityQuad*den
dat$Nc_noClones = dat$Nc*dat$clonalDiversity
dat$NeNc = dat$Ne_SPAGEDI/dat$Nc_noClones
dat$NcFlower = dat$FlowerPlantDensityQuad*den
dat$NcFlower_noClones = dat$NcFlower*dat$clonalDiversity
dat$NeNcFlower = dat$Ne_SPAGEDI/dat$NcFlower_noClones


#print(plot(Nb~Ne_adults_NeEstimator ,dat))
#text(x=dat$Ne_adults_NeEstimator,y=dat$Nb,dat$SvD,cex=.7,pos=4)
#segments(0,0,500,500,lty = "dashed")
#m = lm(Nb~Ne_adults_NeEstimator ,dat)
#abline(m,col="red")
#print(summary(m))

print(bwplot(Ne_SPAGEDI~SvD,dat))
print(wilcox.test(x = dat$Ne_SPAGEDI[dat$SvD=="Shallow"],y = dat$Ne_SPAGEDI[dat$SvD=="Deep"],paired = F))

print(bwplot(NeNc~SvD,dat))
print(wilcox.test(x = dat$NeNc[dat$SvD=="Shallow"],y = dat$NeNc[dat$SvD=="Deep"],paired = F))

print(bwplot(dat$NeNcFlower~SvD,dat))
print(wilcox.test(x = dat$NeNcFlower[dat$SvD=="Shallow"],y = dat$NeNcFlower[dat$SvD=="Deep"],paired = F))


pdf("output/Ne_FieldDensity.pdf",width=5,height=4)
par(mar=c(3,5,3,1))
plot(y=dat$Ne_SPAGEDI, x = as.numeric(factor(dat$SvD)),xlim = c(0.5,2.5),pch=20,cex=2,xlab="",xaxt="n",ylab="Ne - SPAGeDi")
mtext(at=1:2,c("Deep","Shallow"),side=1,line=1.5,cex=1.5)
mtext(side=3,"Wilcox p = 0.8857",line=1.5)

plot(y=dat$NeNc*10^3, x = as.numeric(factor(dat$SvD)),xlim = c(0.5,2.5),pch=20,cex=2,xlab="",xaxt="n",ylab="Nb/Nc (x10^-3)")
mtext(at=1:2,c("Deep","Shallow"),side=1,line=1.5,cex=1.5)
mtext(side=3,"Plant density; Wilcox p = 0.05714",line=1.5)

plot(y=signif(dat$NeNcFlower,3), x = as.numeric(factor(dat$SvD)),xlim = c(0.5,2.5),pch=20,cex=2,xlab="",xaxt="n",ylab="Nb/Nc")
mtext(at=1:2,c("Deep","Shallow"),side=1,line=1.5,cex=1.5)
mtext(side=3,"Flowering plant density; Wilcox p = 0.02857",line=1.5)
dev.off()