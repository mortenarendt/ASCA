# plot kernelsmoother example
rm(list = ls())
library(ggplot2)
obj <- R.matlab::readMat('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA_modelscat/kernelsmoother_example.mat')
X <- data.frame(u = obj$u, obj$X + obj$E, id = 'observations')
Xfun <- data.frame(u = obj$u2, obj$X2,id = 'functional')
Xhat <- data.frame(u = obj$u, obj$Xh,id = 'estimates')

g1 <- bind_rows(X,Xhat) %>% 
  gather(resp,y,X1:X3) %>% 
  ggplot(data = ., aes(u,y, color = resp, shape = id, alpha = id)) + 
  geom_point() + 
  geom_line(data = Xfun %>% gather(resp,y,X1:X3), alpha = 1) + 
  scale_color_brewer(palette = 'Set1') + 
  theme_bw() + 
  ylab('Response') + 
  scale_alpha_manual(values =  c(0.5,1)) + 
  theme(legend.position = 'none',panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text = element_blank())
ggsave('~/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA_modelscat/kernelsmoother_example.pdf',g1, width = 4, height = 4)

