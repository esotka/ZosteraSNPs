### clones on map and calculation of clonal frequency
rm(list=ls())
### meta of adults####
inds.adult <- readLines("data/268adults99seeds_clean")
env <- data.frame(ind=substr(inds.adult,1,7),
                  adult.seed=substr(inds.adult,7,7),
                  site=substr(inds.adult,1,1),
                  SvD=tolower(substr(inds.adult,2,2)),
                  grid=substr(inds.adult,3,3),
                  ind2=substr(inds.adult,5,6))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "cur"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "nil"
env$site.nice[env$site.nice=="W"] <- "wes"
env$site.nice[env$site.nice=="L"] <- "lyn"
env$site.nice <- factor(env$site.nice)
env$quad.name <- paste(env$site.nice,".",env$SvD,env$grid,"_",env$ind2,sep="")
env$site.nice.depth <- paste(env$site.nice,env$SvD,sep="_")
env$site.nice.depth <- factor(env$site.nice.depth)
env$site.nice.depth <- factor(env$site.nice.depth,levels(env$site.nice.depth)[c(1:4,7,8,5,6)])

# Where are the clones? 
# here are geographic positions
clones = read.csv("output/Clone_analysis.csv")
latlon <- read.csv("output/IndivLatLon.csv")
tmp <- ifelse(latlon$ind<10,paste("0",latlon$ind,sep=""),latlon$ind)
latlon$ind <- tmp
latlon$quad.name <- paste(latlon$grid,latlon$ind,sep = "_")
clones$a.quad.name = env$quad.name[match(clones$indA,env$ind)]
clones$b.quad.name = env$quad.name[match(clones$indB,env$ind)]
clones$a.lon = latlon$lon[match(clones$a.quad.name,latlon$quad.name)]
clones$a.lat = latlon$lat[match(clones$a.quad.name,latlon$quad.name)]
clones$b.lon = latlon$lon[match(clones$b.quad.name,latlon$quad.name)]
clones$b.lat = latlon$lat[match(clones$b.quad.name,latlon$quad.name)]

site.nice.depth2 <- c("Curlew_D","Curlew_S","Lynch_D","Lynch_S","West_D","West_S","Niles_D","Niles_S")

pdf("output/ClonesOnMap_Both.pdf",height=6,width=12)
clones2 = clones[clones$adult.seed=="AA" & clones$method=="both",]
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1), mar=c(2,1,0,0))

site.nice.depth2 <- c("Curlew_D","Curlew_S","Lynch_D","Lynch_S","West_D","West_S","Niles_D","Niles_S")

for (i in 1:8)
{
  tmpQuadName <- env$quad.name[env$site.nice.depth==levels(env$site.nice.depth)[i]]
  tmpLat <- latlon$lat[match(tmpQuadName,latlon$quad.name)]; tmpLat = tmpLat[complete.cases(tmpLat)]
  tmpLon <- latlon$lon[match(tmpQuadName,latlon$quad.name)]; tmpLon = tmpLon[complete.cases(tmpLon)]
  plot(tmpLon,tmpLat,ylim=c(min(tmpLat),min(tmpLat)+0.0004),
       xlim=c(min(tmpLon),min(tmpLon)+0.00055))
  mtext(site.nice.depth2[i],line=-2)
  ### colored are clones
  segments(x0=clones2$a.lon,y0=clones2$a.lat,x1=clones2$b.lon,y1=clones2$b.lat,
      col="red",lwd=2)
}

dev.off()  


## how far apart are clones, by site X depth combination
library(geosphere) # geographic distance distm()
m <- c()
for(ii in 1:dim(clones2)[1])
{
m <- c(m,distm(x=clones2[ii,c("a.lon","a.lat")],y=clones2[ii,c("b.lon","b.lat")]))
}
clones2$m <- m

#> summary(clones2$m)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.991   1.998   2.820   2.970   3.982   8.248 

clones2$depth = paste(substr(clones2$indA,2,2),substr(clones2$indB,2,2))
(anova(lm(m~depth,clones2)))
#          Df  Sum Sq Mean Sq F value Pr(>F)
#depth      1   1.877  1.8767  1.0089 0.3193
#Residuals 58 107.884  1.8601

### ## frequency of unique genotypes
allind = env[env$adult.seed=="A",] # 268
(table(allind$site.nice.depth))
clones2$site.nice.depth = allind$site.nice.depth[match(clones2$indA,allind$ind)]
(table(clones2$site.nice.depth))
dat <- 1-table(clones2$site.nice.depth)/table(allind$site.nice.depth)
(dat2 <- data.frame(CD=dat,SvD=rep(c("Deep","Shallow"),4)))
print(anova(lm(CD.Freq~SvD,data=dat2)))
# clonal diversity (1-prop clones)
#  CD.Var1   CD.Freq     SvD
#1   cur_d 0.9069767    Deep
#2   cur_s 0.3684211 Shallow
#3   lyn_d 0.7000000    Deep
#4   lyn_s 0.7500000 Shallow
#5   wes_d 0.9473684    Deep
#6   wes_s 0.8928571 Shallow
#7   nil_d 0.8965517    Deep
#8   nil_s 0.9259259 Shallow
#> print(anova(lm(CD.Freq~SvD,data=dat2)))
#Analysis of Variance Table
#Response: CD.Freq
#          Df   Sum Sq  Mean Sq F value Pr(>F)
#SvD        1 0.032985 0.032985  0.8504  0.392
#Residuals  6 0.232721 0.038787
