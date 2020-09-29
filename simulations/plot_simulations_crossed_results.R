library(tidyverse)
rm(list = ls())
dev.off()
setwd("~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations")
x <- R.matlab::readMat('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/BIGsimulation_crossedPerm.mat')
DOE <- x$DOE %>% data.frame()
pv <- x$pv %>% data.frame()
colnames(pv) <- c('lm','lmm')
colnames(DOE) <- c('k','iter') 
sum(!is.na(pv[,1]))

X <- cbind(DOE,pv)

xx <- X %>% 
  filter(!is.na(lm)) %>% 
  gather(mdl,pv,lm:lmm) %>% 
  group_by(k,mdl) %>% 
  dplyr::summarise(mn = mean(pv), md = median(pv), s = sd(pv))

g1 <- ggplot(data = xx, aes(k,mn,group = mdl,color = mdl))  +
  geom_point() + geom_line() +  
  scale_color_brewer(palette = 'Set1') + theme_bw() + 
  theme(legend.position = 'bottom',
        axis.line = element_line(colour = "black"))
g1

g1
