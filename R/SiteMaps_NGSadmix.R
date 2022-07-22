# site maps with NGSadmix 
rm(list=ls())
library(sf) # geographic distance distm(); destPoint()
library(plotrix) # floating.pie
meta <- as.data.frame(matrix(data=c("cur.s", -70.91575, 42.42009,
                                    "lyn.s", -70.85797, 42.54488,
                                    "wes.s", -70.80600, 42.55922,
                                    "nil.s", -70.65534, 42.59665,
                                    "boston",-71.0589,42.3601),ncol=3,byrow=T))
colnames(meta) <- c("loc","lon","lat")
meta$lon <- as.numeric(meta$lon); meta$lat <- as.numeric(meta$lat)

shp <- st_read("data/NSDE41661/")
shp.df <- as.data.frame(st_coordinates(shp))
tmp <- duplicated(shp.df[,1:2]) # which ones are duplicates
shp.df <- shp.df[tmp==FALSE,]
pdf("output/SiteMaps_NGSadmix.pdf")
plot(shp.df$X,shp.df$Y,col="darkgrey",cex=.1,xlim=c(-71.1,-70.6),ylim=c(42.3,42.65),xlab="",ylab="")

#plot(1,1,xlim=c(-71.1,-70.6),ylim=c(42.3,42.65),xlab="",ylab="")
points(meta$lon,meta$lat,col="black",pch=c(20,20,20,20,21),cex=3)
# scale bar
leftPoint <- c(-70.75,42.35)
rightPoint <- destPoint(leftPoint,90,5000)
segments(x0=leftPoint[1],y0=leftPoint[2],x1=rightPoint[1],y1=rightPoint[2],col="black",lwd=2)
text("5 km",x=(leftPoint[1]+rightPoint[1])/2,y=leftPoint[2]+0.01,cex=.8)
text(x=meta$lon,y=meta$lat+c(-0.02,0.02,0.03,0.02,-0.02),c("Curlew Beach","Lynch Park","West Beach","Niles Beach","Boston MA"))
#text("Atlantic Ocean",x=(leftPoint[1]+rightPoint[1])/2,y=42.4,cex=1.5,font=3)
arrows(x0=-71.05,y0=42.55,x1=-71.05,y1=42.60,length=0.2,col="black",lwd=3)
text("N",x=-71.05,y=42.60,pos=3,cex=1.5)


## admix k=3 (highest delta). No clones
a <- read.delim("data/ngsadmix-noClones/k03run1.qopt",header=F,sep=" ")
kcol <- c("red","darkgreen","black")
inds <- readLines("data/ind245adults_noClones_clean")
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
env$site.nice <- factor(env$site.nice,levels=levels(env$site.nice)[c(1,2,4,3)])


xbar <- aggregate(a[,1:3],by=list(env$site.nice),mean)
# pie position
piex <- c(-70.91037,-70.84924,-70.78069,-70.67046)
piey <- c(42.34091,42.41549,42.46572,42.47789)
segments(meta$lon[1:4],meta$lat[1:4],piex,piey)
for (j in 1:4)
{floating.pie(xpos=piex[j],ypos=piey[j],x=c(t(xbar[j,-1])),col=kcol,radius=0.02)}

# seeds
s <- read.delim("data/NGSadmix-seeds/k03run1.qopt",header=F,sep=" ")
inds <- readLines("data/ind99_seeds")
env <- data.frame(ind=substr(inds,1,7),
                  adult.seed=substr(inds,7,7),
                  site=substr(inds,1,1),
                  SvD=substr(inds,2,2))
## change D to C (Dorothy ==> Curlew Beach)
env$site.nice <- as.character(env$site)
env$site.nice[env$site.nice=="D"] <- "Curlew"  #Dorothy Cove
env$site.nice[env$site.nice=="N"] <- "Niles"
env$site.nice[env$site.nice=="W"] <- "West"
env$site.nice[env$site.nice=="L"] <- "Lynch"
env$site.nice <- factor(env$site.nice)
env$site.nice <- factor(env$site.nice,levels=levels(env$site.nice)[c(1,2,4,3)])
xbar <- aggregate(s[,1:3],by=list(env$site.nice),mean)
for (j in 1:4)
{floating.pie(xpos=piex[j]+0.05,ypos=piey[j],x=c(t(xbar[j,-1])),col=kcol,radius=0.02)}
text(x=piex,y=piey,"A",col="white")
text(x=piex+0.05,y=piey,"S",col="white")
dev.off()
