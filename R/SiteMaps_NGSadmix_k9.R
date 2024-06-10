# site maps with NGSadmix 
rm(list=ls())
library(sf) # geographic distance distm(); st_read
library(geosphere) #destPoint()
library(plotrix) # floating.pie
meta <- as.data.frame(matrix(data=c("Curlew", -70.91575, 42.42009,-70.91037,42.34,
                                    "Lynch", -70.85797, 42.54488,-70.75,42.395,
                                    "West", -70.80600, 42.55922,-70.7,42.4705,
                                    "Niles", -70.65534, 42.59665,-70.67046,42.54,
                                    "Boston",-71.0589,42.3601,NA,NA),ncol=5,byrow=T))
# pie position
colnames(meta) <- c("loc","lon","lat","piex","piey")
row.names(meta) = meta$loc; meta = meta[,-1]
for (i in 1:4){meta[,i] = as.numeric(meta[,i])}
#meta$lon <- as.numeric(meta$lon); meta$lat <- as.numeric(meta$lat)

shp <- st_read("data/NSDE41661/")
shp.df <- as.data.frame(st_coordinates(shp))
tmp <- duplicated(shp.df[,1:2]) # which ones are duplicates
shp.df <- shp.df[tmp==FALSE,]
pdf("output/SiteMaps_NGSadmix_k9.pdf")
plot(shp.df$X,shp.df$Y,col="darkgrey",cex=.1,xlim=c(-71.1,-70.6),ylim=c(42.3,42.65),xlab="",ylab="")

points(meta$lon,meta$lat,col="black",pch=c(20,20,20,20,21),cex=3)
text(x=meta$lon,y=meta$lat+0.02,cex=1.6,rownames(meta))
text("Atlantic Ocean",x=-70.68,y=42.32,cex=1.5,font=3)
meta = meta[1:4,] ## remove boston
# segments to pies
segments(meta$lon[1:4],meta$lat[1:4],
       x1=c(meta$piex[1],meta$piex[2:3]-0.06,meta$piex[4]),
       y1=meta$piey)

# scale bar
leftPoint <- c(-71.08,42.5)
rightPoint <- destPoint(leftPoint,90,5000)
segments(x0=leftPoint[1],y0=leftPoint[2],x1=rightPoint[1],y1=rightPoint[2],col="black",lwd=2)
text("5 km",x=(leftPoint[1]+rightPoint[1])/2,y=leftPoint[2]+0.01,cex=.8)
arrows(x0=-71.05,y0=42.55,x1=-71.05,y1=42.60,length=0.2,col="black",lwd=3)
text("N",x=-71.05,y=42.60,pos=3,cex=1.5)


## admix k=5 (highest delta). No clones
a <- read.delim("data/ngsAdmixFiles/k09run1.qopt",header=F,sep=" ")
k=9
#kcol <- c("red","darkgreen","black","gainsboro","yellow")[c(4,3,5,1,2)]
#kcol <- c("red","darkgreen","black","gainsboro","yellow","deepskyblue","brown")[c(6,1,7,4,3,5,2)]
kcol <- c("black", #1
          "red",#2
          "darkgreen",#3
          "gainsboro",#4
          "yellow",#5
          "deepskyblue",#6
          "brown",#7
          "dodgerblue4",#8
          "darkorchid")[c(7,6,2,1,8,9,5,3,4)]#9
inds <- readLines("data/214adults99seeds_noClones_clean")
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
print(levels(env$site.nice))
#[1] "Curlew" "Lynch"  "West"   "Niles"
xbar <- aggregate(a[,1:k],by=list(env$site.nice,env$adult.seed,env$SvD),mean)
xbar$pos = meta[match(xbar$Group.1,rownames(meta)),]
# adultshallow
tmp = xbar[xbar$Group.2=="A" & xbar$Group.3=="S",]
text("AS",x=tmp$pos$piex,y=tmp$pos$piey-0.03,col="black")
for (j in 1:4)
{floating.pie(xpos=tmp$pos$piex[j],ypos=tmp$pos$piey[j],x=c(t(tmp[j,4:(k+3)])),col=kcol,radius=0.025)}
    
# adultdeep
tmp = xbar[xbar$Group.2=="A" & xbar$Group.3=="D",]
text("AD",x=tmp$pos$piex-0.06,y=tmp$pos$piey-0.03,col="black")
for (j in 1:4)
{floating.pie(xpos=tmp$pos$piex[j]-0.06,ypos=tmp$pos$piey[j],x=c(t(tmp[j,4:(k+3)])),col=kcol,radius=0.025)}

xbar <- aggregate(a[,1:k],by=list(env$site.nice,env$adult.seed),mean)
xbar$pos = meta[match(xbar$Group.1,rownames(meta)),]
# all seeds
tmp = xbar[xbar$Group.2=="S",]
text("Seed",x=tmp$pos$piex+0.06,y=tmp$pos$piey-0.03,col="black")
for (j in 1:4)
{floating.pie(xpos=tmp$pos$piex[j]+0.06,ypos=tmp$pos$piey[j],x=c(t(tmp[j,3:(k+2)])),col=kcol,radius=0.025)}


#text(x=c(meta$piex-0.06,meta$piex,meta$piex+0.06),y=meta$piey,c("DA","SA","S"),col="white")
#text(x=piex+0.05,y=piey+0.027,"S",col="black")


dev.off()
