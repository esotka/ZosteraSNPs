### seed phenotype
rm(list=ls())
library(ggplot2)
library(gridExtra)
dat = read.csv("data/seedPhenotypes.csv")
dat$mass = as.numeric(dat$mass) # normally-distributed

### number of seed per core
out = aggregate(dat$seedTF,by=list(dat$site,dat$depth,dat$quad,dat$core),sum)
colnames(out) = c(colnames(dat)[1:4],"NumSeeds")

## binomial model - seed or no seed
m1 = glm(out$NumSeeds==0 ~ out$site * out$depth, family = binomial)
print(anova(m1,test="Chisq"))
e1 = resid(m1,type="pearson")
dispersion = sum(e1^2) / m1$df.resid
print(dispersion)
#0.8414634 very little overdispersion

# there's a significant interaction between site and depth;
# stats by site
for (i in 1:4)
{
    tmp = out[out$site==unique(out$site)[i],]
    print(unique(tmp$site))
    m1 = glm(tmp$NumSeeds==0 ~ tmp$depth, family = binomial)
    print(anova(m1,test="Chisq"))

}
# DC p = 0.002844
# LP p = 0.2642
# NB p = 0.07812
# WB p = 3.046e-07

print(aggregate(out$NumSeeds,by=list(out$site,out$depth),mean))

#  Group.1 Group.2          x
#1      DC      DP 0.55555556
#2      LP      DP 0.08333333
#3      NB      DP 0.16666667
#4      WB      DP 0.00000000
#5      DC      SH 8.55555556
#6      LP      SH 0.41666667
#7      NB      SH 0.91666667
#8      WB      SH 5.75000000

meta <- data.frame(
  site = c("DC","LP","WB","NB"),
  site.nice = c("Curlew","Lynch","West","Niles"))
out$site.nice <- meta$site.nice[match(out$site,meta$site)]
out$site.nice <- factor(out$site.nice)
out$site.nice <- factor(out$site.nice,levels=levels(out$site.nice)[c(1,2,4,3)])
out$siteDepth <- paste(out$site.nice,out$depth,sep="_")


f1 <- ggplot(out, aes(x=site.nice, y=NumSeeds, fill=depth)) +
  geom_boxplot() + 
  ylab("Number of seeds per core") + xlab("") + ylim(0,15) + 
  theme_classic(base_size=20) + theme(legend.position=c(0.2,0.9)) +
  scale_fill_grey(start=.75,end=1,name="") +
  annotate(geom = "text", x = c(1,2,3,4), y = c(7.5,5,15,5), label = c("0.002","ns","<0.001","ns"),size=5)

### mass per seed

seed = dat[dat$seedTF,]
seed$site.nice <- meta$site.nice[match(seed$site,meta$site)]
seed$site.nice <- factor(seed$site.nice)
seed$site.nice <- factor(seed$site.nice,levels=levels(seed$site.nice)[c(1,2,4,3)])
seed$siteDepth <- paste(seed$site.nice,seed$depth,sep="_")


m2 = lm(mass~site,seed)
print(anova(m2))

#Response: mass
#           Df     Sum Sq    Mean Sq F value    Pr(>F)    
#site        3 8.5387e-05 2.8462e-05  12.028 5.469e-07 ***
#Residuals 128 3.0290e-04 2.3664e-06

print(aggregate(seed$mass,by=list(seed$site),mean,na.rm=T))

#  Group.1           x
#1      DC 0.008343750
#2      LP 0.006000000
#3      NB 0.005791667
#4      WB 0.007343396

print(TukeyHSD(aov(m2)))

# DC A
# WB B
# LP ABC
# NB C

#               diff           lwr           upr     p adj
#LP-DC -0.0023437500 -0.0047092365  2.173653e-05 0.0531255
#NB-DC -0.0025520833 -0.0038117623 -1.292404e-03 0.0000033
#WB-DC -0.0010003538 -0.0017440553 -2.566522e-04 0.0035268
#NB-LP -0.0002083333 -0.0027931400  2.376473e-03 0.9967192
#WB-LP  0.0013433962 -0.0010330563  3.719849e-03 0.4576724
#WB-NB  0.0015517296  0.0002715768  2.831882e-03 0.0106406


f2 <- ggplot(seed, aes(x=site.nice, y=mass)) +
  geom_boxplot() + 
  ylab("grams per seed") + xlab("") + #+ ylim(0,15) + 
  theme_classic(base_size=20) + theme(legend.position=c(0.2,0.9)) +
  scale_fill_grey(start=.75,end=1,name="") +
  annotate(geom = "text", x = c(1,2,3,4), y = 0.012, label = c("A","ABC","B","C"),size=5)

png('output/seedPhenotypes.png',width=12,height=7,units="in",res=700)
grid.arrange(f1,f2,nrow=1,ncol=2)
dev.off()
