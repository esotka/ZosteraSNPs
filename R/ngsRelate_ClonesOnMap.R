### do ngsRelate values (r>0.8) show up as clones? 

rm(list=ls())
### meta of adults####
inds.adult <- readLines("data/ind294adults_clean")
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

### ngsRelate data
a <- c("DD",
       "DS",
       "LD",
       "LS",
       "ND",
       "NS",
       "WD",
       "WS")

rab <- c()
for (i in 1:length(a))
{
  tmp <- read.delim(file=paste("data/ngsRelate/",a[i],"_newres",sep=""))[,c("a","b","rab")]
  tmpind <- readLines(paste("data/ngsRelate/",a[i],sep=""))
  tmp.quad.name <- env$quad.name[match(tmpind,env$ind)]
  tmp.key <- data.frame(num=0:(length(tmpind)-1),tmp.quad.name)
  tmp$a <- tmp.key$tmp.quad.name[match(tmp$a,tmp.key$num)]
  tmp$b <- tmp.key$tmp.quad.name[match(tmp$b,tmp.key$num)]
  rab <- rbind(rab,tmp)
}

rab$site <- substr(rab$a,1,3)
rab$depth <- substr(rab$a,5,5)
rab$site_depth <- paste(rab$site,rab$depth,sep="_")
library(lattice)
print(histogram(~rab,data=rab,col="grey",breaks=40))
print(histogram(~rab,data=rab,col="grey",breaks=40,ylim=c(0,5)))

# Where are the clones? 
# here are geographic positions
latlon <- read.csv("output/IndivLatLon.csv")
tmp <- ifelse(latlon$ind<10,paste("0",latlon$ind,sep=""),latlon$ind)
latlon$ind <- tmp
latlon$quad.name <- paste(latlon$grid,latlon$ind,sep = "_")

threshold <- 0.8 ### the threshold for relatedness. 
rab.sub <- rab[rab$rab>=threshold,]
rab.sub$a.lon <- latlon$lon[match(rab.sub$a,latlon$quad.name)]
rab.sub$a.lat <- latlon$lat[match(rab.sub$a,latlon$quad.name)]
rab.sub$b.lon <- latlon$lon[match(rab.sub$b,latlon$quad.name)]
rab.sub$b.lat <- latlon$lat[match(rab.sub$b,latlon$quad.name)]

pdf("output/ngsRelate_ClonesOnMap=threshold=0.8.pdf",height=6,width=12)
par(mfcol=c(2,4),mai=c(0.3,0.3,0.1,0.1), mar=c(2,1,0,0))

site.nice.depth2 <- c("Curlew_D","Curlew_S","Lynch_D","Lynch_S","West_D","West_S","Niles_D","Niles_S")

for (i in 1:8)
{
  tmpQuadName <- env$quad.name[env$site.nice.depth==levels(env$site.nice.depth)[i]]
  tmpLat <- latlon$lat[match(tmpQuadName,latlon$quad.name)]
  tmpLon <- latlon$lon[match(tmpQuadName,latlon$quad.name)]
  plot(tmpLon,tmpLat,ylim=c(min(tmpLat),min(tmpLat)+0.0004),
       xlim=c(min(tmpLon),min(tmpLon)+0.00055))
  mtext(site.nice.depth2[i],line=-2)
  ### colored are clones
  segments(x0=rab.sub$a.lon,y0=rab.sub$a.lat,x1=rab.sub$b.lon,y1=rab.sub$b.lat,col="red",lwd=2)
}
  
dev.off()  

## how far apart are clones, by site X depth combination
library(geosphere) # geographic distance distm()
m <- c()
for(ii in 1:dim(rab.sub)[1])
{
m <- c(m,distm(x=rab.sub[ii,c("a.lon","a.lat")],y=rab.sub[ii,c("b.lon","b.lat")]))
}
rab.sub$m <- m
#Group.1   x.Min. x.1st Qu. x.Median   x.Mean x.3rd Qu.   x.Max.
#1   cur_d 1.991228  1.991228 1.995614 2.201947  2.206333 2.825331
#2   cur_s 1.991228  2.000000 2.825331 3.343264  4.117802 8.248334
#3   lyn_d 1.991228  1.996519 2.006105 2.725374  2.828427 6.325939
#4   lyn_s 1.991228  1.991228 1.999141 2.367865  2.821178 3.982447
#5   nil_d 1.998283  1.999571 2.991224 2.990794  3.982447 3.982447
#6   nil_s 1.998283  3.080197 4.162111 4.162111  5.244025 6.325939
#7   wes_s 1.993895  1.998283 2.006105 2.541243  2.826105 4.004386

  
### New list of unique genotypes 
# manually #
a.alt <- env$ind[match(rab.sub$a,env$quad.name)]
b.alt <- env$ind[match(rab.sub$b,env$quad.name)]

out <- list()
out[[1]] <- c(a.alt[1],b.alt[1])
for (k in 2:length(a.alt))
{
  a.query <- sapply(out, function(y) a.alt[k] %in% y)
  b.query <- sapply(out, function(y) b.alt[k] %in% y)
    if(sum(a.query) > 0)
    {
      tmp <- unique(c(unlist(out[a.query]),b.alt[k]))
      out[[(1:length(out))[a.query]]] <- tmp
    }
  else
  {
    if(sum(b.query)>0)
      {
      tmp <- unique(c(unlist(out[b.query]),a.alt[k]))
      out[[(1:length(out))[b.query]]] <- tmp
    }
    else
    {
      out[[k]] <- c(a.alt[k],b.alt[k])
    }
  }
}

### get rid of zeros
out <- out[lengths(out)>0]
sink("output/ListOfClones.txt")
print(out)
sink()

## frequency of unique genotypes
clones <- unlist(out)
nonClonalgenets <- env$ind[!(env$ind%in%clones)] # pull out all clones
genets <- unlist(lapply(out,"[[",1)) # keep the first one
uniqueGenets <- c(nonClonalgenets,genets)
print(table(substr(uniqueGenets,1,2))) # unique genotype numbers
print(table(substr(env$ind,1,2))) # all individuals sequenced
print(dat <- table(substr(uniqueGenets,1,2))/table(substr(env$ind,1,2)))
dat2 <- data.frame(CD=dat,SvD=rep(c("Deep","Shallow"),4))
print(anova(lm(CD.Freq~SvD,data=dat2)))


### made a list of adults with clones pulled. ***Do not change***

#clones.to.remove <- c()
#clones.to.keep <- c()
#for (m in 1:length(out))
#{
#  tmp <- out[[m]] # keep a randomly chosen clone
#  tmp.to.remove <- tmp[(sample(1:length(tmp),size = length(tmp)-1))]
#  clones.to.remove <- c(clones.to.remove,tmp.to.remove)
#  tmp.to.keep <- tmp[!(tmp%in%tmp.to.remove)]
#  clones.to.keep <- c(clones.to.keep,tmp.to.keep)
#}

### length(clones.to.remove) = 49
### length(clones.to.keep) = 32
### noClones <- env$ind[!(env$ind%in%clones.to.remove)]
### writeLines(noClones,"data/ind245adults_noClones_clean")

