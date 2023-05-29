## Maps of sites
rm(list=ls())
library(geosphere) # geographic distance distm(); destPoint()
library(scales) # alpha()
pdf("output/SiteMaps.pdf",width=10,height=10)
par(mfrow=c(2,2))
# Curlew Beach = Dorothy Cove
cur.d3 <- c(-70.91754,42.41904) #measured
cur.d2 <- destPoint(cur.d3,78.09183,18) #  heading; 18 meters away
cur.d1 <- destPoint(cur.d2,78.09183,18) # heading; 18 meters away
#cur.d2 <- c(-70.91722,42.41909) #measured
# bearing(cur.d3,cur.d2) yields 78.09183
#cur.d1 <- c(-70.91690,42.41914) ### this was not linear, so I used the slope for the next pt. 
#lat <- cur.d2[1]+abs(cur.d3[1]-cur.d2[1])
#lon <- cur.d2[2]+(abs(cur.d3[2])-abs(cur.d2[2]))
cur.s1 <- c(-70.91553,42.42009)
cur.s2 <- destPoint(cur.s1,270,18) # 270º heading; 18 meters away
cur.s3 <- destPoint(cur.s2,270,18) # 270º heading; 18 meters away


cur <- data.frame(lon=c(cur.d3[1],cur.d2[1],cur.d1[1],
                        cur.s1[1],cur.s2[1],cur.s3[1]),
                lat=c(cur.d3[2],cur.d2[2],cur.d1[2],
                      cur.s1[2],cur.s2[2],cur.s3[2]),
                row.names=c("cur.d3","cur.d2","cur.d1",
                            "cur.s1","cur.s2","cur.s3"))
plot(cur$lon,cur$lat,xlab="longitude",ylab="latitude",type="n",main="Curlew")
text(x=cur$lon,y=cur$lat,substr(row.names(cur),5,6))

### Niles Beach

nil.d1 <- c(-70.65647,42.59595)
nil.d2 <- destPoint(nil.d1,270,12)
nil.d3 <- destPoint(nil.d1,270,36)
nil.s1 <- c(-70.65523,42.59654)
nil.s2 <- c(-70.65549,42.59668)
nil.s3 <- c(-70.65531,42.59674)

nil <- data.frame(lon=c(nil.d3[1],nil.d2[1],nil.d1[1],
                        nil.s1[1],nil.s2[1],nil.s3[1]),
                  lat=c(nil.d3[2],nil.d2[2],nil.d1[2],
                        nil.s1[2],nil.s2[2],nil.s3[2]),
                  row.names=c("nil.d3","nil.d2","nil.d1",
                              "nil.s1","nil.s2","nil.s3"))
plot(nil$lon,nil$lat,xlab="longitude",ylab="latitude",type="n",main="Niles")
text(x=nil$lon,y=nil$lat,substr(row.names(nil),5,6))

### West beach

wes.d1 <- c(-70.80425,42.55788)
wes.d2 <- destPoint(wes.d1,-86.48149,18)
wes.d3 <- destPoint(wes.d2,-86.48149,18)
#wes.d2 <- c(-70.80456,42.55785)
#wes.d3 <- c(-70.80469,42.55790)
# bearing(wes.d1,wes.d3) ==> -86.48149
wes.s1 <- c(-70.80578,42.55921)
wes.s2 <- destPoint(wes.s1,-86.48149,18)
wes.s3 <- destPoint(wes.s2,-86.48149,18)

wes <- data.frame(lon=c(wes.d3[1],wes.d2[1],wes.d1[1],
                        wes.s1[1],wes.s2[1],wes.s3[1]),
                  lat=c(wes.d3[2],wes.d2[2],wes.d1[2],
                        wes.s1[2],wes.s2[2],wes.s3[2]),
                  row.names=c("wes.d3","wes.d2","wes.d1",
                              "wes.s1","wes.s2","wes.s3"))
plot(wes$lon,wes$lat,xlab="longitude",ylab="latitude",type="n",main="West Beach")
text(x=wes$lon,y=wes$lat,substr(row.names(wes),5,6))

### Lynch

lyn.d1 <- c(-70.85712,42.54358) 
lyn.d2 <- destPoint(lyn.d1,62.01922,18)
lyn.d3 <- destPoint(lyn.d2,62.01922,18)
#lyn.d2 <- c(-70.85701,42.54362) 
#lyn.d3 <- c(-70.85684,42.54369) 
# bearing(lyn.d1,lyn.d3) ==> 62.01922
lyn.s1 <- c(-70.85816,42.54480) 
lyn.s2 <- destPoint(lyn.s1,62.01922,18)
lyn.s3 <- destPoint(lyn.s2,62.01922,18)

lyn <- data.frame(lon=c(lyn.d3[1],lyn.d2[1],lyn.d1[1],
                        lyn.s1[1],lyn.s2[1],lyn.s3[1]),
                  lat=c(lyn.d3[2],lyn.d2[2],lyn.d1[2],
                        lyn.s1[2],lyn.s2[2],lyn.s3[2]),
                  row.names=c("lyn.d3","lyn.d2","lyn.d1",
                              "lyn.s1","lyn.s2","lyn.s3"))
plot(lyn$lon,lyn$lat,xlab="longitude",ylab="latitude",type="n",main="Lynch Beach")
text(x=lyn$lon,y=lyn$lat,substr(row.names(lyn),5,6))


all <- rbind(cur,nil,wes,lyn)
write.csv(all,"output/SiteQuadratLatLon.csv")
library(maps)
library(mapdata)
library(maptools) 
#library(scales)  
map("worldHires","USA",xlim=c(-71,-70),ylim=c(42,43),col="gainsboro",fill=TRUE)
points(all$lon,all$lat,pch=20,col="red")#,xlab="longitude",ylab="latitude")
hist(distm(x=all,fun=distGeo)*10^-3,xlab="km",main="km distance between grids")

### within grid distances
### 2 m between each point.
library(sp)
g <- GridTopology(c(0,0), c(2,2), c(5,3))
g.coord <- coordinates(g)
g.coord.std <- data.frame(s1=g.coord[,1]-4,s2=g.coord[,2]-2) # standardize for center point being c(0,0)
g.coord.std[,3]  <- c(1,12,11,10,9,2,13,14,15,8,3:7)
colnames(g.coord.std) <- c("s1","s2","ind")
plot(g.coord.std[,-3],type="n")
text(g.coord.std[,1],g.coord.std[,2],g.coord.std[,3])
print(g.coord.std)

## coordinates of points 2m away from center of grid
g.coord.std$metersDist <- sqrt(g.coord.std$s1^2+g.coord.std$s2^2) # pythagorean therom (a2 + b2 = c2)
# orientation of points away from center of grid

### Assume all sites are on the same plane (90º) 
#  will this matter? ==> Curlew Deep axis is 78.09183
g.coord.std$lon <- NA; g.coord.std$lat <- NA
g2 <- g.coord.std[order(g.coord.std$ind),] # sort by number

### get bearings: use these grid numbers as lat and lon
bearing.out <- c()
for (i in 1:15)
  {bearing.out <- c(bearing.out,bearing(g2[14,c("s1","s2")],g2[i,c("s1","s2")]))}

g2$bearing <- bearing.out

### generate sample lat and lon
tmp <- all[rownames(all)=="cur.d3",]
g2[14,c("lon","lat")] <- tmp ### center of the grid
for (i in c(1:13,15))
{g2[i,c("lon","lat")] <- destPoint(g2[14,c("lon","lat")],
                                   g2$bearing[i],g2[i,"metersDist"])}

plot(g2$lon,g2$lat,type="n",main="ind samples at Curlew Deep 3")
text(g2$lon,g2$lat,g2$ind)
#hist(distm(g2[,c("lon","lat")]),breaks=20,main="m between samples at Curlew Deep grid 3",xlab="m")

################################################
### make a list of all sample lat / lon ###
################################################
all.ind <- c()
gridIDs <- rownames(all)
for (j in 1:length(gridIDs))
{

tmp <- all[rownames(all)==gridIDs[j],]
g2[14,c("lon","lat")] <- tmp ### center of the grid
for (i in c(1:13,15))
{g2[i,c("lon","lat")] <- destPoint(g2[14,c("lon","lat")],
                                   g2$bearing[i],g2[i,"metersDist"])}

all.ind <- rbind(all.ind,data.frame(grid=gridIDs[j],
                                    pop=substr(gridIDs[j],1,3),
                                    ind=g2$ind,
                                    lon=g2$lon,
                                    lat=g2$lat))
}

write.csv(all.ind,"output/IndivLatLon.csv")

### print them
pops.unique <- unique(all.ind$pop)
for (k in 1:4)
{
tmp <- all.ind[all.ind$pop==pops.unique[k],]
plot(tmp$lon,tmp$lat,pch=20,cex=.5,xlab="",ylab="",main=pops.unique[k])
}

################################################
### make a distance matrix in meters
################################################

tmp <- ifelse(all.ind$ind<10,paste("0",all.ind$ind,sep=""),all.ind$ind)
all.ind$ind <- tmp

all.dist <- distm(all.ind[,c("lon","lat")])
par(mfrow=c(1,1))
hist(all.dist,xlab="meters",breaks=50)
dev.off()

rownames(all.dist) <- paste(all.ind$grid,all.ind$ind,sep = "_")
colnames(all.dist) <- paste(all.ind$grid,all.ind$ind,sep = "_")

all.dist2 <- as.data.frame(all.dist)
write.csv(all.dist2,"output/GeographicDistanceIndiv.csv",quote=F)

#### distance between depths ####

#lon.xbar lat.xbar
#cur.d -70.91733 42.41907
#cur.s -70.91575 42.42009
#lyn.d -70.85693 42.54366
#lyn.s -70.85797 42.54488
#nil.d -70.65666 42.59595
#nil.s -70.65534 42.59665
#wes.d -70.80447 42.55789
#wes.s -70.80600 42.55922

#habitat distances: 
#  Curlew: 170m
#  Lynch:  160m
#  Niles:  133m
#  West:   193m

#site distances: Cur, Lyn, Nil, wes (in kilometers)

 #        [,1]      [,2]     [,3]      [,4]
#[1,]  0.00000 14.677733 29.04944 17.945792
#[2,] 14.67773  0.000000 17.52077  4.572644
#[3,] 29.04944 17.520769  0.00000 12.949062
#[4,] 17.94579  4.572644 12.94906  0.000000

################################################
### make a pretty map with a coastline ###
################################################
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
cols.unique <- c(alpha(c("blue","purple","red","black"),.7),"black")
points(meta$lon,meta$lat,col=cols.unique,pch=c(20,20,20,20,21),cex=3)
# scale bar
leftPoint <- c(-70.75,42.35)
rightPoint <- destPoint(leftPoint,90,5000)
segments(x0=leftPoint[1],y0=leftPoint[2],x1=rightPoint[1],y1=rightPoint[2],col="black",lwd=2)
text("5 km",x=(leftPoint[1]+rightPoint[1])/2,y=leftPoint[2]+0.01,cex=.8)
text(x=meta$lon,y=meta$lat+c(-0.02,0.02,-0.02,0.02,-0.02),
     c("Curlew Beach","Lynch Park","Niles Beach","West Beach","Boston MA"),
     col=cols.unique)
text("Atlantic Ocean",x=(leftPoint[1]+rightPoint[1])/2,y=42.4,cex=1.5,font=3)
arrows(x0=-71.05,y0=42.55,x1=-71.05,y1=42.60,length=0.2,col="black",lwd=3)
text("N",x=-71.05,y=42.60,pos=3,cex=1.5)
dev.off()



