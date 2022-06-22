### Fis with adults vs seed

rm(list=ls())
library(ggplot2)
dat <- read.csv("output/BasicStats_adultsNoClones.csv")
dat$site_depth <- paste(dat$site,dat$depth)

anova(lm(Fis~depth+site+adult.seed,dat))
#Response: Fis
#Df    Sum Sq   Mean Sq F value Pr(>F)  
#depth       1 0.0000059 0.0000059  0.0032 0.9568  
#site        3 0.0062400 0.0020800  1.1064 0.4083  
#adult.seed  1 0.0158734 0.0158734  8.4431 0.0228 *
#  Residuals   7 0.0131603 0.0018800                 

pdf("output/Fis~adult.seed.pdf",width=5,height=6)
f <- ggplot(dat, aes(x=adult.seed, y=Fis, fill="darkgray")) +
  geom_boxplot(outlier.size = NULL) +
  ylab("Fis") + xlab("") +
  theme_classic(base_size=20) + theme(legend.position="none") +
  scale_fill_grey(start=.75,end=1)
print(f)
dev.off()


