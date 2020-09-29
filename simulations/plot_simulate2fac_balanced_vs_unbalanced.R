setwd("~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA_modelscat")
rm(list = ls())
library(ggplot2)
library(tidyverse)
xx <- R.matlab::readMat('simulate2fac_balanced_vs_unbalanced_v2.mat')
df <- data.frame(q = xx$qq  %>% as.vector(),
                 pv1 = xx$pv1 %>% as.vector(),
                 pv2 = xx$pv2%>% as.vector())

# calculate mean, median and variance and compare with theoretical
res <- df %>% 
gather(meth,pv,pv1,pv2) %>% 
  group_by(q,meth) %>% 
  dplyr::summarise(mn = mean(pv), 
            vr = var(pv))

res %>% 
  mutate(biasmn = mn - 0.5, 
         biasvr = vr - 1/12) %>% 
  gather(metric,bias,biasmn, biasvr) %>% 
  ggplot(data = ., aes(q,bias,color = meth, group = meth)) + 
  geom_line() + facet_wrap(~metric,nrow = 1, scales = 'free') + 
  theme_bw() + geom_hline(yintercept = 0)




xx2 <- R.matlab::readMat('testing_conditionalpermutation.mat')
df <- data.frame(xx2$RES)
colnames(df) <- c('N','corr obs','mean bias','median bias')
g1 <- df %>% 
  mutate(corr = round(0.5*`corr obs`,1)*2) %>% 
  mutate(N = factor(N)) %>% 
  group_by(N,corr) %>% 
  summarise(mean = mean(`mean bias`),
            median = mean(`median bias`)) %>% 
  gather(meth,bias,mean,median) %>% 
  ggplot(data = ., aes(corr,bias, color = N, linetype = meth )) + 
  geom_line() +
  geom_point() + 
  scale_color_brewer(palette = 'Set2') + 
  theme_bw() + theme(legend.position = 'bottom')  

pdf('Simulations_condiperm.pdf')
g1
dev.off()

df <- data.frame(q = xx$qq  %>% as.vector(),
                 pv1 = xx$pv1 %>% as.vector(),
                 pv2 = xx$pv2%>% as.vector())


xx <- R.matlab::readMat('testing_propmatching.mat')
df <- data.frame(xx$RES)
colnames(df) <- unlist(xx$lbs)
df$orth <- 'Orthogonalized'

xx2 <- R.matlab::readMat('testing_propmatching_woOrth.mat')
df2 <- data.frame(xx2$RES)
colnames(df2) <- unlist(xx2$lbs)
df2$orth <- 'Not orthogonalized'

head(df)
bind_rows(df,df2) %>% 
  mutate(propmatch = ifelse(propmatch==1,'Mathcing','No Matching'), 
         prop_orth = paste(propmatch,orth,sep = ' - ')) %>% 
  group_by(Q,orth,prop_orth) %>%
  dplyr::summarise(mean = mean(pv),
            median = median(pv)) %>% 
  gather(meth,pv,mean,median) %>% 
  ggplot(data = ., aes(Q,pv, color = prop_orth, group = prop_orth)) + 
  geom_line() +
  geom_point() + 
  facet_wrap(~meth) + ylim(c(0,1)) + 
  # scale_color_brewer(palette = 'Set2') + 
  theme_bw() + theme(legend.position = 'bottom')  
