### ngsRelate results - seeds

a <- c("DS_seeds",
       "WS_seeds",
       "NS_seeds")

## within grid vs between grid
library(ggplot2)
rab <- c()
for (i in 1:length(a))
{
  tmp <- read.delim(file=paste("data/ngsRelate/",a[i],"_newres",sep=""))[,c("a","b","rab")]
  tmpind <- readLines(paste("data/ngsRelate/",a[i],sep=""))
#  tmp.quad.name <- substr(tmpind,3,3)
  tmp.key <- data.frame(num=0:(length(tmpind)-1),tmpind)
  tmp$a <- tmp.key$tmpind[match(tmp$a,tmp.key$num)]
  tmp$b <- tmp.key$tmpind[match(tmp$b,tmp.key$num)]
  rab <- rbind(rab,tmp)
}
rab$gridA <- substr(rab$a,3,3)
rab$gridB <- substr(rab$b,3,3)
rab$gridAB <- paste(rab$gridA,rab$gridB)
rab$site <- substr(rab$a,1,1)
rab$site[rab$site=="D"] <- "C"
rab$site <- factor(rab$site); rab$site <- factor(rab$site,levels=levels(rab$site)[c(1,3,2)])
rab$WithinBetweenGrid <- ifelse(rab$gridAB%in%c("1 1","2 2","3 3"),"withinGrid","betweenGrid")  
  

print(table(rab$gridAB,rab$site))

#      D   N   W
#1 1   1   1  78
#1 2  16  14 221
#1 3  50   2 208
#2 2  28  21 136
#2 3 200   7 272
#3 3 300   0 120

# remove groups that do not have > 5 seeds

rab.sub <- rab[!(rab$gridAB=="1 1" & rab$site=="D"),]
rab.sub <- rab.sub[!(rab.sub$gridAB=="1 1" & rab.sub$site=="N"),]
rab.sub <- rab.sub[!(rab.sub$gridAB=="1 3" & rab.sub$site=="N"),]
rab.sub <- rab.sub[!(rab.sub$gridAB=="2 3" & rab.sub$site=="N"),]

library(lattice)
print(bwplot(rab~gridAB | site,data=rab.sub,col="grey"))
print(bwplot(rab~WithinBetweenGrid | site,data=rab.sub,col="grey"))

mod <- glm(rab~WithinBetweenGrid*site,data=rab.sub,family="poisson")
library(AER)
dispersiontest(mod)
anova(mod)

#Terms added sequentially (first to last)


#Df Deviance Resid. Df Resid. Dev
#NULL                                    1663     173.51
#WithinBetweenGrid       1   46.454      1662     127.05
#site                    2    2.251      1660     124.80
#WithinBetweenGrid:site  2    0.210      1658     124.59

library(pgirmess)
#print(PermTest(mod,B=1000))
#PermTest.glm(obj = mod, B = 1000)

#Based on 1000 replicates
#Simulated p-value:
#  p.value
#WithinBetweenCore        0.000
#site                     0.000
#WithinBetweenCore:site   0.450


#### core
library(readxl)
meta <- read_xlsx("data/Seed ID to core.xlsx")
meta$grid_core <- paste(meta$grid,meta$core_within_grid,sep="_")
rab$grid_coreA <- meta$grid_core[match(substr(rab$a,1,6),substr(meta$RAD_name,1,6))]
rab$grid_coreB <- meta$grid_core[match(substr(rab$b,1,6),substr(meta$RAD_name,1,6))]
rab$grid_coreAB <- paste(rab$grid_coreA,rab$grid_coreB,sep=" ")

tmp <-  c(paste(paste(1:4,LETTERS[1],sep="_"),paste(1:4,LETTERS[1],sep="_"),sep=" "),
          paste(paste(1:4,LETTERS[2],sep="_"),paste(1:4,LETTERS[2],sep="_"),sep=" "),
          paste(paste(1:4,LETTERS[3],sep="_"),paste(1:4,LETTERS[3],sep="_"),sep=" "),
          paste(paste(1:4,LETTERS[4],sep="_"),paste(1:4,LETTERS[4],sep="_"),sep=" "))

rab$WithinBetweenCore <- ifelse(rab$grid_coreAB%in%tmp,"withinCore","betweenCore")  

rab$GridCoreCombo <- paste(rab$WithinBetweenGrid,rab$WithinBetweenCore)

rab$GridCoreCombo[rab$GridCoreCombo=="betweenGrid betweenCore"] <- "betweenGrid"
rab$GridCoreCombo[rab$GridCoreCombo=="withinGrid betweenCore"] <- "betweenCore"
rab$GridCoreCombo[rab$GridCoreCombo=="withinGrid withinCore"] <- "withinCore"
rab$GridCoreCombo <- factor(rab$GridCoreCombo)
rab$GridCoreCombo <- factor(rab$GridCoreCombo,levels=levels(rab$GridCoreCombo)[c(2,1,3)])

print(table(rab$grid_coreAB,rab$site))
#print(bwplot(rab~grid_coreAB | site,data=rab,col="grey"))
print(bwplot(rab~GridCoreCombo| site,data=rab,col="grey"))
pdf("output/Relatedness~Core-ngsRelate_seeds.pdf")
f <- ggplot(rab, aes(x=site, y=rab, fill=GridCoreCombo)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("rab") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="right",legend.title=element_blank()) +
  scale_fill_grey(start=.75,end=1)
print(f)
mod <- glm(rab~GridCoreCombo*site,data=rab,family="poisson")
#Df Deviance Resid. Df Resid. Dev
#NULL                                1674     174.91
#GridCoreCombo       2   48.387      1672     126.53
#site                2    2.249      1670     124.28
#GridCoreCombo:site  4    0.583      1666     123.69
print(anova(mod))
#print(PermTest(mod,B=1000))
#Based on 1000 replicates
#Simulated p-value:
#  p.value
#GridCoreCombo         0.000
#site                  0.000
#GridCoreCombo:site    0.334

f <- ggplot(rab[rab$GridCoreCombo%in%c("withinCore","betweenCore"),], aes(x=site, y=rab, fill=GridCoreCombo)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("rab") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="right",legend.title=element_blank()) +
  scale_fill_grey(start=.75,end=1)
print(f)
dev.off()
print(bwplot(rab~GridCoreCombo| site,data=rab[rab$GridCoreCombo%in%c("withinCore","betweenCore"),],col="grey"))
mod2 <- glm(rab~GridCoreCombo*site,data=rab[rab$GridCoreCombo%in%c("withinCore","betweenCore"),],family="poisson")
print(anova(mod2))
#Df Deviance Resid. Df Resid. Dev
#NULL                                 684     76.441
#GridCoreCombo       1  1.22599       683     75.215
#site                2  1.93821       681     73.276
#GridCoreCombo:site  2  0.29059       679     72.986
#print(PermTest(mod2,B=1000))
#Based on 1000 replicates
#Simulated p-value:
#  p.value
#GridCoreCombo         0.000
#site                  0.000
#GridCoreCombo:site    0.294

