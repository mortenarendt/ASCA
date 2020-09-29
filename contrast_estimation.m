%%
% 
%  Simulated example on the use of hierachical contrast estimation using
%  ASCA

%% Simulate data
p = 20; % number of variables
k = 6; % Factor1
rep = 10; % Factor2
D = fullfact([rep k]);

Dtrue = D;
Dtrue(ismember(Dtrue(:,1),[1 2 3]),1) = 1; % merge the three first levels into one. 

%% MODEL
M = [1 0; 0 1]; % additive model

[X DM Xm] = simulateData(Dtrue,M,[0.1 0.1 .2],p,[]);

%% ASCA settings
plotit=0;
opt = ASCAcat('options');
opt.nperm = 1000;
opt.showtable = 0;
opt.rndfac = [];

%% Compute ASAC model
res1= ASCAcat(X,D,M, opt);
res1.ANOVAtab

%% Compute hieracical contrasts
rescont = ASCAcontrast(res1,1);
plotASCAcontrast(rescont); 
