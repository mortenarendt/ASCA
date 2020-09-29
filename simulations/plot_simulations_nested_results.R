library(tidyverse)
library(scales)
rm(list = ls())
setwd("~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations")
x <- R.matlab::readMat('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/BIGsimulation_nestedPerm2.mat')
DOE <- x$DOE %>% data.frame()
pv <- x$pv %>% data.frame()
S2N <- x$S2N %>% as.vector() %>% data.frame()
B2N <- x$background2N %>% as.vector() %>% data.frame()
colnames(S2N) <- 'S2N'
colnames(B2N) <- 'B2N'
colnames(pv) <- c('lm','lmm')
colnames(DOE) <- c('k','iter') 
sum(!is.na(pv[,1]))

X <- cbind(DOE,S2N,B2N,pv) %>% data.frame()
xx <- X %>% 
  filter(!is.na(lm)) %>% 
  gather(mdl,pv,lm:lmm) %>% 
  group_by(k,mdl) %>% 
  dplyr::summarise(mn = mean(pv), md = median(pv), s = sd(pv),q1 = quantile(pv)[2],q3 = quantile(pv)[4], 
                   S2N = mean(S2N), B2N = mean(B2N)) 

w <- 0.3
g1 <- xx %>% 
  filter(k>0) %>% 
  ggplot(data = ., aes(S2N,mn,ymin = q1, ymax =q3,group = mdl,color = mdl))  +
  geom_point(position = position_dodge(width = w)) + 
  # geom_point() + 
  geom_line(position = position_dodge(width = w)) +
  geom_errorbar(width = 0.5,position = position_dodge(width = w)) +
  # geom_point() + geom_line() +  
  geom_point(data = xx[xx$k==0,], aes(x = 5e-7),position = position_dodge(width = w)) + 
  geom_errorbar(data = xx[xx$k==0,], aes(x = 5e-7),position = position_dodge(width = w)) + 
  # scale_x_log10() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_log10() + 
  ylab('p') + 
  scale_color_brewer(palette = 'Set1') + theme_bw() + 
  theme(legend.position = 'none',
        axis.line = element_line(colour = "black"))
g1


pdf('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/sim_nested.pdf', width = 3, height = 3)
g1
dev.off()
