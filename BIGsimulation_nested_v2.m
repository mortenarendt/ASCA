close all hidden; clear; clc
cd '/Users/mortenarendtrasmussen/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/'

% test LMM for nested designs, where the nested factor is null, but the
% coarser increase in level (x-axis)

%%% settings
nperm = 50; % number of permutations in the ASCA model
Slinspan = [0.01 0.1];
REP = [1 2 5 10];
niter = 50; % for each setting, run niter models

DOE = fullfact([length(Slinspan) length(REP) niter]);
DOE(:,1) = Slinspan(DOE(:,1));
DOE(:,2) = REP(DOE(:,2));

% set seed
rng(321)

p = 20;
lbD = {'ID','condition', 'trt'};
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
% MODEL
M = [1 0 0 ; 0 1 0 ; 0 0 1];
s = [0.1 0 0.1 .2];
plotit=0;
opt = ASCAcat('options');
opt.nperm = nperm;
opt.showtable = 0;

%%% RECORD
tic
pv = nan(size(DOE,1),2);
for j=1:size(DOE,1)
    j
    rep = DOE(j,2); % repeated measures
    D1 = fullfact([rep k]);
    D2 = (D1(:,2)<16)+1;
    D = [D1(:,2) D2 D1(:,1)]; % ID, condition, trt
    
    s(1) = DOE(j,1);
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
    
    % save results
    if round(j / 100)==(j/100)
        save('./simulations/BIGsimulation_nestedPerm_v3.mat','DOE','pv','S2N','background2N')
    end
end

save('./simulations/BIGsimulation_nestedPerm_v3.mat','DOE','pv','S2N','background2N')
toc

%
%for i=1:length(Slinspan)
%    p = pv(DOE(:,1)==Slinspan(i),:);
%    PP(i,:) = mean(p,1);
%end
%pv
%plot(Slinspan,PP);    shg

