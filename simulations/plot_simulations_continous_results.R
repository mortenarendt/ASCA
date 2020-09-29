library(cowplot)
library(tidyverse)
library(scales) # to access break formatting functions
rm(list = ls())

x <- R.matlab::readMat('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/BIGsimulation_combined.mat')
DOE <- x$DOE %>% data.frame()
pv <- x$pvall %>% data.frame()
EXPVAR <- x$EXPVARall %>% data.frame()
colnames(EXPVAR)[-1] <- colnames(pv) <- c('lin','KS','KSslc')
colnames(pv) <- paste('pv_',colnames(pv),sep = '')
colnames(EXPVAR)[1] <- 'truth'
colnames(DOE) <- c('k','iter') 

X <- cbind(DOE,EXPVAR,pv)

x <- R.matlab::readMat('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/BIGsimulation_continous3.mat')
DOE <- x$DOE %>% data.frame()
colnames(DOE) <- c('k','iter') 
EXPVAR <- x$EXPVAR %>% data.frame()
RELFIT <- x$RELFIT %>% data.frame()
CORR <- x$CORR %>% data.frame()
S2N <- x$S2N %>% t() %>% data.frame()
colnames(CORR) <- colnames(RELFIT) <- colnames(EXPVAR)[-1] <- c('lin','KS','KSslc')
colnames(EXPVAR)[1] <- 'truth'
colnames(S2N) <- 'S2N'
X2 <- cbind(DOE,S2N,EXPVAR)
X3 <- cbind(DOE,S2N,truth = EXPVAR[,1],CORR)

S2Ndf <- X2 %>% 
  group_by(k) %>% 
  dplyr::summarise(S2N = mean(S2N))

xx <- X2 %>% 
  filter(!is.na(lin)) %>%
  filter(k > 0) %>% 
  gather(mdl,EV,lin:KSslc) %>% 
  mutate(dEV = EV / truth) %>% 
  # mutate(dEV = EV) %>%
  group_by(k,mdl) %>% 
  dplyr::summarise(mn = mean(dEV), md = median(dEV), q1 = quantile(dEV)[2],q3 = quantile(dEV)[4]) %>% 
  left_join(S2Ndf, by = 'k')

w <- 0.3
g1 <- ggplot(data = xx, aes(S2N,mn,ymin = q1, ymax =q3,
                            group = mdl,color = mdl))  +
  geom_point(position = position_dodge(width = w)) + 
  geom_line(position = position_dodge(width = w)) + 
  geom_errorbar(width = 0.5,position = position_dodge(width = w)) + 
  # scale_x_log10() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_log10() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_brewer(palette = 'Set1') + theme_bw() + 
  ylab('SSQ(marginal) / SSQ(truth)') + 
  theme(legend.position = 'none',
        axis.line = element_line(colour = "black"))

g1
xx3 <- X3 %>% 
  filter(!is.na(lin)) %>%
  filter(k > 0) %>% 
  gather(mdl,r,lin:KSslc) %>% 
  group_by(k,mdl) %>% 
  # dplyr::summarise(mn = mean(r), md = median(r), s = sd(r)) %>% 
  dplyr::summarise(mn = mean(r), md = median(r), q1 = quantile(r)[2],q3 = quantile(r)[4]) %>% 
  left_join(S2Ndf, by = 'k')

g3 <- ggplot(data = xx3, aes(S2N,mn,ymin = q1, ymax =q3,
                             group = mdl,color = mdl))  +
  geom_point(position = position_dodge(width = w)) + 
  geom_line(position = position_dodge(width = w)) + 
  geom_errorbar(width = 0.5,position = position_dodge(width = w)) + 
                             
  # scale_x_log10() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_brewer(palette = 'Set1') + theme_bw() + 
  ylab('r') + 
  theme(legend.position = 'none',
        axis.line = element_line(colour = "black"))





xx2 <- X %>% 
  filter(!is.na(lin)) %>% 
  gather(mdl,p,pv_lin:pv_KSslc) %>% 
  group_by(k,mdl) %>% 
  # dplyr::summarise(mn = mean(p), md = median(p), s = sd(p)) %>% 
  dplyr::summarise(mn = mean(p), md = median(p), q1 = quantile(p)[2],q3 = quantile(p)[4]) %>% 
  left_join(S2Ndf, by = 'k')

brks <- 10^seq(-8,-1,2)
lbs <- paste()
g2 <- xx2 %>% filter(S2N>0) %>% 
  ggplot(data = ., aes(S2N,mn,ymin = q1, ymax =q3,
                             group = mdl,color = mdl))  +
  geom_point(position = position_dodge(width = w)) + 
  # geom_point() + 
  geom_line(position = position_dodge(width = w)) +
  geom_errorbar(width = 0.5,position = position_dodge(width = w)) +
  
  geom_point(data = xx2[xx2$S2N==0,], aes(x = 1e-7),position = position_dodge(width = w)) + 
  geom_errorbar(data = xx2[xx2$S2N==0,], aes(x = 1e-7),position = position_dodge(width = w)) + 
  
  coord_cartesian(ylim = c(0,1)) +
  # scale_x_log10(breaks = brks) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab('p') + 
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.line = element_line(colour = "black"))

pdf('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/sim_continous.pdf', width = 8, height = 3)
plot_grid(g1, g3, g2, labels = c('A', 'B','C'), label_size = 12, nrow = 1)
dev.off()
