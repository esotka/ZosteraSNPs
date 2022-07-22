### Fst with distance - no Clone adults
rm(list=ls())
library(geosphere) # geographic distance distm()

# phist (from amova_noClones.R)
dat <- read.csv("output/amova.method1-calls_noClones-pairwisePhiSt.csv")

# get pairwise distance
pos <- data.frame(matrix(c(
  "Curlew_D","-70.91733","42.41907",
  "Curlew_S","-70.91575","42.42009",
  "Lynch_D","-70.85693","42.54366",
  "Lynch_S","-70.85797","42.54488",
  "West_D","-70.80447","42.55789",
  "West_S","-70.80600","42.55922",
  "Niles_D","-70.65666","42.59595",
  "Niles_S","-70.65534","42.59665"
  ),nrow=8,ncol=3,byrow = T))
colnames(pos) <- c("site", "lon","lat")
geo <- distm(pos[,c("lon","lat")])
colnames(geo) <- pos$site; rownames(geo) <- pos$site

dat$km <- c()
for (i in 1:dim(dat)[1])
{dat$km[i] <- geo[match(dat$pop1[i],rownames(geo)),match(dat$pop2[i],colnames(geo))]/1000}



h1 <- unlist(lapply(strsplit(dat$pop1,"_"),"[[",2))
h2 <- unlist(lapply(strsplit(dat$pop2,"_"),"[[",2))
dat$hab <- factor(paste(h1,h2,sep=""))

cols <- c("red","black","black","blue")
pdf("output/Fst~Distance.pdf",width=7,height=5)
plot(dat$km,dat$phiST,xlab="Geographic distance (km)",ylab="phiSt",pch=19,col=cols[dat$hab])
abline(lm(dat$phiST~(dat$km)))
arrows(1,.107,.5,.085)
text(x=0,y=0.11,pos=4,"Between habitat")
dev.off()

cor.test(dat$phiST,dat$km,method="kendall")

dat$hab[dat$hab=="SD"] <- "DS"
dat$hab <- factor(dat$hab)
anova(lm(phiST~hab,dat)) # no significant difference between DD, DS, SS

#Analysis of Variance Table

#Response: phiST
#Df    Sum Sq    Mean Sq F value Pr(>F)
#hab        2 0.0021507 0.0010753  1.0599 0.3616
#Residuals 25 0.0253647 0.0010146
anova(lm(phiST~hab,dat[!dat$hab=="DS",])) # no significant difference between DD or SS

#Analysis of Variance Table

#Response: phiST
#Df    Sum Sq    Mean Sq F value Pr(>F)
#hab        1 0.0013659 0.00136589  1.8154 0.2076
#Residuals 10 0.0075239 0.00075239              

## we don't have much power to detect differences in the slope of these (low sample size)

