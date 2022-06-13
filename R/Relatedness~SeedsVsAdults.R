## ngsRelate seeds vs adults 

## within grid vs between grid
## seeds
library(lattice)
library(readxl)
rm(list=ls())
a <- c("DS_seeds",
       "WS_seeds",
       "NS_seeds")

seed <- c()
for (i in 1:length(a))
{
  tmp <- read.delim(file=paste("data/ngsRelate/",a[i],"_newres",sep=""))[,c("a","b","rab")]
  tmpind <- readLines(paste("data/ngsRelate/",a[i],sep=""))
  #  tmp.quad.name <- substr(tmpind,3,3)
  tmp.key <- data.frame(num=0:(length(tmpind)-1),tmpind)
  tmp$a <- tmp.key$tmpind[match(tmp$a,tmp.key$num)]
  tmp$b <- tmp.key$tmpind[match(tmp$b,tmp.key$num)]
  seed <- rbind(seed,tmp)
}
seed$gridA <- substr(seed$a,3,3)
seed$gridB <- substr(seed$b,3,3)
seed$gridAB <- paste(seed$gridA,seed$gridB)
seed$site <- substr(seed$a,1,1)
seed$WithinBetweenGrid <- ifelse(seed$gridAB%in%c("1 1","2 2","3 3"),"withinGrid","betweenGrid")  


print(table(seed$gridAB,seed$site))
#      D   N   W
#1 1   1   1  78
#1 2  16  14 221
#1 3  50   2 208
#2 2  28  21 136
#2 3 200   7 272
#3 3 300   0 120

# remove groups that do not have > 5 seeds

seed.sub <- seed[!(seed$gridAB=="1 1" & seed$site=="D"),]
seed.sub <- seed.sub[!(seed.sub$gridAB=="1 1" & seed.sub$site=="N"),]
seed.sub <- seed.sub[!(seed.sub$gridAB=="1 3" & seed.sub$site=="N"),]
seed.sub <- seed.sub[!(seed.sub$gridAB=="2 3" & seed.sub$site=="N"),]

# adults - clones removed
#inds.adult <- readLines("data/ind245adults_noClones_clean")
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

a <- c("DD",
       "DS",
      # "LD", # no data in seeds
      # "LS", # no data in seeds
       "ND",
       "NS",
       "WD",
       "WS")

adult <- c()
for (i in 1:length(a))
{
  tmp <- read.delim(file=paste("data/ngsRelate/",a[i],"_newres",sep=""))[,c("a","b","rab")]
  tmpind <- readLines(paste("data/ngsRelate/",a[i],sep=""))
  tmp.quad.name <- env$quad.name[match(tmpind,env$ind)]
  tmp.key <- data.frame(num=0:(length(tmpind)-1),tmp.quad.name)
  tmp$a <- tmp.key$tmp.quad.name[match(tmp$a,tmp.key$num)]
  tmp$b <- tmp.key$tmp.quad.name[match(tmp$b,tmp.key$num)]
  adult <- rbind(adult,tmp)
}

adult$gridA <- substr(adult$a,5,6)
adult$gridB <- substr(adult$b,5,6)
adult$gridAB <- paste(adult$gridA,adult$gridB)
adult$site <- substr(adult$a,1,3)
adult$site[adult$site=="cur"] <- "D"
adult$site[adult$site=="nil"] <- "N"
adult$site[adult$site=="wes"] <- "W"
adult$WithinBetweenGrid <- ifelse(adult$gridAB%in%c("d1 d1","d2 d2","d3 d3","s1 s1","s2 s2","s3 s3"),"withinGrid","betweenGrid")  

### remove clones

noClones <- readLines("data/ind245adults_noClones_clean")
env$noClones <- env$ind%in%noClones
adult$noClonesA <- env$noClones[match(adult$a,env$quad.name)]
adult$noClonesB <- env$noClones[match(adult$b,env$quad.name)]
adult_noClones <- adult[(adult$noClonesA+adult$noClonesB)==2,] # both are unique


### combine
alldat <- rbind(seed,adult_noClones[,1:8])
alldat$AvS <- c(rep("seed",dim(seed)[1]),rep("adult",dim(adult_noClones)[1]))

print(bwplot(rab~ AvS | WithinBetweenGrid ,data=alldat,col="grey"))
print(bwplot(rab~ AvS | site+WithinBetweenGrid ,data=alldat,col="grey"))

#mod <- lm(log(rab+0.000001)~site+WithinBetweenGrid*AvS,data=alldat)#,family="poisson")
mod <- glm(rab~site+WithinBetweenGrid*AvS,data=alldat,family="poisson")

library(AER)
dispersiontest(mod)
print(anova(mod))


#library(pgirmess)
#print(PermTest(mod,B=1000))

