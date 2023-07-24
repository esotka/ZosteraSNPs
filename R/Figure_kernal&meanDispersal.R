### visualize dispersal kernals and mean dispersal distance estimates
### with help from Pinsky et al. 2016
library(scales)
rm(list=ls())
dat = read.table("output/sigma_estimates.txt",header=T)
# Assuming a Gaussian kernal (following Siegel et al. 2003 MEPS equation 6)

meanDispersal = dat$sigma_point*sqrt(2/pi); print(meanDispersal)
meanDispersal_L95 = dat$sigmaL95*sqrt(2/pi); print(meanDispersal_L95)
meanDispersal_U95 = dat$sigmaU95*sqrt(2/pi); print(meanDispersal_U95)

# using 200m linear distance
#meanDispersal = dat$sigma_point*sqrt(2/pi); print(meanDispersal)
#[1] 0.3306963 0.2846050 0.2497856 0.2339052
#> meanDispersal_L95 = dat$sigmaL95*sqrt(2/pi); print(meanDispersal_L95)
#[1] 0.2824511 0.2425569 0.2130352 0.1994711
#> meanDispersal_U95 = dat$sigmaU95*sqrt(2/pi); print(meanDispersal_U95)
#[1] 0.4180915 0.3582502 0.3143665 0.2944194

# using 500m linear distance
#meanDispersal 
#[1] 0.1653482 0.1423025 0.1248928 0.1169526
#meanDispersal_L95
#[1] 0.14122557 0.12127845 0.10691653 0.09973557
#meanDispersal_U95
#[1] 0.2090458 0.1795240 0.1571833 0.1476086 

# Gaussian kernal (following Siegel et al. 2003 MEPS equation 5, where xd = 0). 

cols.to.use = grey.colors(4)
out = c()
pdf("output/kernal&meanDispersal_Figure.pdf",width=8,height=5)
par(fig = c(0,1,0,1))
plot(x=0:2,y=seq(0,1,.5),type="n",xlab="Distance (km)",ylab="Dispersal strength (normalized)",xlim=c(0,5))
for (i in 1:4)
{
tmp = dat[i,c("sigma_point","sigmaL95","sigmaU95")]
x = seq(0,10,0.01)
dispersalKernal = function(x,tmp){(1/(tmp*sqrt(2*pi)))*exp((-x^2)/(2*tmp))}
px = dispersalKernal(x,as.numeric(tmp[1])); px = px/max(px) # mean
px_L95 = dispersalKernal(x,as.numeric(tmp[2])); px_L95 = px_L95/max(px_L95) # L95
px_U95 = dispersalKernal(x,as.numeric(tmp[3])); px_U95 = px_U95/max(px_U95) # L95
points(x,px,type="l",col=cols.to.use[i])
polygon(c(x,sort(x,decreasing = T)),c(px_L95,sort(px_U95)),col=alpha(cols.to.use[i],.2))
}

legend(x=1,y=.8,dat$MAF,fill=cols.to.use,bty = "n")
text(1.35,0.9,"MAF",cex=1.3,pos=1)

par(fig = c(.55,.95,.3,.95), new = T)
plot(meanDispersal~c(4:1),type="p",xaxt="n",
     ylim=c(0,0.25),xlim=c(0.5,4.5),pch=20,cex.lab=0.8,xlab="",ylab="")
arrows(4:1,meanDispersal_L95,4:1,meanDispersal_U95,angle = 90,code = 3,length = 0.1)
mtext(at=4:1,side = 1,line=-1,dat$MAF)
mtext("Mean dispersal (km)",side=2,line=2.5)
dev.off()


