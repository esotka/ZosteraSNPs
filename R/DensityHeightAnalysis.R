#####################################################
## Analysis of shoot morphology in permanent quadrats 
## 25 cm x 25 cm quadrats set at every other node 
## spatially paired with our genetic sampling 
#####################################################
library(lattice)
library(lmerTest)
library(ggplot2)
library(car) # Anova
library(gridExtra) # gridArrange
rm(list=ls())
dat.all <- read.csv("data/perm.quad.density.2019.csv")
# parse to include just the June data from the grids; the August data was collected when we collected flowers for the outplant experiment later that year; not spatially paired with SNPs data, but same rough location 

#################################################
#### Total density (=vegetative + flowering) ####
#################################################

dat <- dat.all[dat.all$month=="June",]
dat$depth <- ifelse(dat$depth=="DP","Deep","Shallow")
dat$biomass <- as.numeric(dat$biomass)
meta <- data.frame(
  site = c("DC","LP","WB","NB"),
  site.nice = c("Curlew","Lynch","West","Niles"))
dat$site.nice <- meta$site.nice[match(dat$site,meta$site)]
dat$site.nice <- factor(dat$site.nice)
dat$site.nice <- factor(dat$site.nice,levels=levels(dat$site.nice)[c(1,2,4,3)])
dat$siteDepth <- paste(dat$site.nice,dat$depth,sep="_")

# linear model - residusals are not 
print(histogram(~total.density | site.nice+depth,data = dat,col="grey",breaks=20))
m1 <- lm(total.density~factor(site)*factor(depth),dat)
print(Anova(m1))

htmodel <- lmer(total.density ~ site * depth + (1|perm.quadrat), data = dat, REML = FALSE)
print(Anova(htmodel))

print(aggregate(dat$total.density,by=list(dat$site.nice,dat$depth),mean))


f1 <- ggplot(dat, aes(x=site.nice, y=total.density, fill=depth)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Plants / 0.0625 m2") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position=c(0.2,0.9)) +
  scale_fill_grey(start=.5,end=1,name="")
print(f1)

#######################################################
#### Aboveground Biomass (=vegetative + flowering) ####
#######################################################

print(histogram(~biomass | site+depth,data = dat,col="grey",breaks=20))
print(Anova(lm(biomass~factor(site)*factor(depth),dat)))

print(aggregate(dat$biomass,by=list(dat$site.nice,dat$depth),mean,na.rm=T))
htmodel <- lmer(biomass ~ site * depth + (1|perm.quadrat), data = dat, REML = FALSE)
print(Anova(htmodel))

f2 <- ggplot(dat[complete.cases(dat$biomass),], aes(x=site.nice, y=biomass, fill=depth)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Aboveground drymass (g) / 0.0625 m2") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="none") +  
  scale_fill_grey(start=.75,end=1)
print(f2)

#######################################################
#### Proportion flowering ####
#######################################################

## zero truncated - used poisson
print(histogram(~f.density | site+depth,data = dat,col="grey",breaks=20))
print(bwplot(f.density~factor(siteDepth),data=dat))
m1 <- glm(f.density~depth+site,dat,family=poisson)
print(Anova(m1))

E.model1<-resid(m1)
Fit.model1<-fitted(m1)
plot(x=Fit.model1, y=E.model1, xlab="fitted values", ylab="vegq_residuals")

f3 <- ggplot(dat, aes(x=site.nice, y=f.density, fill=depth)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Flowering plants / 0.0625 m2") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="none") +  
  scale_fill_grey(start=.75,end=1)
print(f3)

## Is there proportionally more Flowering in the deep or shallow? No.
print(histogram(~perc.flowering | site+depth,data = dat,col="grey",breaks=20))
m1 <- glm(perc.flowering~depth*site,dat,family=poisson)
Anova(m1)

fx <- ggplot(dat, aes(x=site.nice, y=perc.flowering, fill=depth)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Proportion of flowering plants") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="none") +  
  scale_fill_grey(start=.75,end=1)
print(fx)



#######################################################################
#### Plant height - up to five lengths measured for each quadrat ######
#######################################################################

dat.all <- read.csv("data/perm.quad.veg.ht.2019.csv")
dat <- dat.all[dat.all$Date=="June",]
dat$depth <- ifelse(dat$depth=="DP","Deep","Shallow")
meta <- data.frame(
  site = c("DC","LP","WB","NB"),
  site.nice = c("Curlew","Lynch","West","Niles"))
dat$site.nice <- meta$site.nice[match(dat$site,meta$site)]
dat$site.nice <- factor(dat$site.nice)
dat$site.nice <- factor(dat$site.nice,levels=levels(dat$site.nice)[c(1,2,4,3)])
dat$siteDepth <- paste(dat$site.nice,dat$depth,sep="_")

print(aggregate(dat$height,by=list(dat$site.nice,dat$depth),mean,na.rm=T))

out <- aggregate(dat$height,by=list(dat$site.nice,dat$depth,dat$unique.quad),mean)
out <- data.frame(out)
colnames(out) <- c("Site","Depth","Quadrat","Height")
# stats
print(histogram(~height | site+depth,data = dat,col="grey",breaks=20))
htmodel <- lmer(height ~ site * depth + (1|unique.quad), data = dat, REML = FALSE)
print(Anova(htmodel))

# posthoc - by site

for (i in 1:4)
{
  tmp <- dat[dat$site==unique(dat$site)[i],]
  htmodel <- lmer(height ~ depth + (1|unique.quad), data = tmp, REML = FALSE)
  print(unique(tmp$site)); print(Anova(htmodel))
}


f4 <- ggplot(out, aes(x=Site, y=Height, fill=Depth)) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Plant height (cm)") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="none") +  
  scale_fill_grey(start=.75,end=1)
print(f4)

png('output/DensityHeight.png',width=12,height=12,units="in",res=700)
grid.arrange(f1,f2,f3,f4,nrow=2,ncol=2)
dev.off()


