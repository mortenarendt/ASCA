clear; clc; close all;
cd '/Users/mortenarendtrasmussen/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA_modelscat'
load BeerData.mat
xs = myscale(lgX.data); 

F = [lgX.class{1,1}' lgX.class{1,2}']; 
Flb = {lgX.classname{1,1};lgX.classname{1,2}}; 
M = [1 0; 0 1]; 
m2 = ASCAcat(xs,F,Flb,100,M);
figure; plotASCAscores(m2,2)
%%
% close all; 
figure
[u s v] = svds(xs,2); 
subplot(2,2,1); scatter(u(:,1),u(:,2),100,lgX.class{1,1},'filled'); title(['PCA (' lgX.classname{1,1} ')'])
text(u(:,1),u(:,2),lgX.label{1})

subplot(2,2,2); scatter(u(:,1),u(:,2),100,lgX.class{1,2},'filled'); title(['PCA (' lgX.classname{1,2} ')'])
text(u(:,1),u(:,2),lgX.label{1})
u = m2.Effects{1}.loads{1}; 
subplot(2,2,3); scatter(u(:,1),u(:,2),100,lgX.class{1,1},'filled'); title(['ASCA (' lgX.classname{1,1} ')'])
u = m2.Effects{2}.loads{1}; 
subplot(2,2,4); scatter(u(:,1),u(:,2),100,lgX.class{1,2},'filled'); title(['ASCA (' lgX.classname{1,2} ')'])
shg; 
% PrintandMerge('PCA_vs_ASCA.pdf');
%%
movefile('PCA_vs_ASCA.pdf','/Users/mortenarendt/Dropbox (Huttenhower Lab)/Backup/MyDocumentsOnC/MATLAB/work/TRAC_chemometrics_in_foodomics/PCA_vs_ASCA.pdf'      )
%%



