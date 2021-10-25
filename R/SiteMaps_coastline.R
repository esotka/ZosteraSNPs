# Coastline Map
# Shapefile taken from https://www.ngs.noaa.gov/CUSP/
rm(list=ls())

library(geosphere)
library(sf)

meta <- as.data.frame(matrix(data=c("cur.s", -70.91575, 42.42009,
           "lyn.s", -70.85797, 42.54488,
           "nil.s", -70.65534, 42.59665,
           "wes.s", -70.80600, 42.55922,
           "boston",-71.0589,42.3601),ncol=3,byrow=T))
colnames(meta) <- c("loc","lon","lat")
meta$lon <- as.numeric(meta$lon); meta$lat <- as.numeric(meta$lat)

shp <- st_read("data/NSDE41661/")
shp.df <- as.data.frame(st_coordinates(shp))
tmp <- duplicated(shp.df[,1:2]) # which ones are duplicates
shp.df <- shp.df[tmp==FALSE,]
pdf("output/SiteMaps_coastline.pdf")
plot(shp.df$X,shp.df$Y,col="darkgrey",cex=.1,xlim=c(-71.1,-70.6),ylim=c(42.3,42.65),xlab="",ylab="")
points(meta$lon,meta$lat,col="black",pch=c(20,20,20,20,21),cex=3)
# scale bar
leftPoint <- c(-70.75,42.35)
rightPoint <- destPoint(leftPoint,90,5000)
segments(x0=leftPoint[1],y0=leftPoint[2],x1=rightPoint[1],y1=rightPoint[2],col="black",lwd=2)
text("5 km",x=(leftPoint[1]+rightPoint[1])/2,y=leftPoint[2]+0.005,cex=.8)
text(x=meta$lon,y=meta$lat+c(-0.02,0.02,-0.02,0.02,-0.02),c("Curlew Beach","Lynch Park","Niles Beach","West Beach","Boston MA"))
text("Atlantic Ocean",x=(leftPoint[1]+rightPoint[1])/2,y=42.4,cex=1.5,font=3)
arrows(x0=-71.05,y0=42.55,x1=-71.05,y1=42.60,length=0.2,col="black",lwd=3)
text("N",x=-71.05,y=42.60,pos=3,cex=1.5)
dev.off()
