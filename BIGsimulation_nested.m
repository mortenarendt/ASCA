close all hidden; clear; clc
cd '/Users/mortenarendtrasmussen/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/'

% test LMM for nested designs, where the nested factor is null, but the
% coarser increase in level (x-axis)

%%% settings
nperm = 50; % number of permutations in the ASCA model
Slinspan = [0 exp(linspace(log(1e-7),log(1e-4),10))];
niter = 50; % for each setting, run niter models

length(Slinspan)*niter;
DOE = fullfact([length(Slinspan) niter]);
DOE(:,1) = Slinspan(DOE(:,1));

% set seed
rng(321)

p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
D = [D1(:,2) D2 D1(:,1)]; % ID, condition, trt
lbD = {'ID','condition', 'trt'};

% make unbalanced
% remove 50% of samples D2==1 and D1 in 1 and 2
% remove 50% of samples D2==2 and D1 in 4 and 5
icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>1;
icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>1;
icout = icout1 | icout2;
D = D(~icout,:);

% MODEL
M = [1 0 0 ; 0 1 0 ; 0 0 1];
s = [0.1 0 0.1 .2];
%s(1) = Slinspan

plotit=0;
opt = ASCAcat('options');
opt.nperm = nperm;
opt.showtable = 0;


%%% RECORD
pv = nan(size(DOE,1),2);
tic
for j=1:size(DOE,1)
    %%
    opt.nperm = 1000; 
    j
    
    s(1) = DOE(j,1);
    s(1) = 1e-12
    % simulate data
    [X DM Xm] = simulateData(D,M,s,p,opt.contvar);
    
    % models 1 - Vanilla (wrong)
    opt.rndfac = [];
    res1= ASCAcat(X,D(:,2:3),M(2:3,2:3), opt);
    
    % models 2 - nested permutation
    opt.rndfac = 1;
    res2 = ASCAcat(X,D,M, opt);
    
    
    %%% export results
    pv(j,:) = [res1.Effects{1}.p res2.Effects{2}.p];
    
    S2N(j) = trace(Xm{1}'*Xm{1}) / trace(Xm{end}'*Xm{end});
    XX = Xm{2} + Xm{3};
    background2N(j) = trace(XX'*XX) / trace(Xm{end}'*Xm{end});
    %%
    % save results
    if round(j / 100)==(j/100)
        save('./simulations/BIGsimulation_nestedPerm2.mat','DOE','pv','S2N','background2N')
    end
end

save('./simulations/BIGsimulation_nestedPerm2.mat','DOE','pv','S2N','background2N')
toc

%
%for i=1:length(Slinspan)
%    p = pv(DOE(:,1)==Slinspan(i),:);
%    PP(i,:) = mean(p,1);
%end
%pv
%plot(Slinspan,PP);    shg

