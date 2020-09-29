close all hidden; clear; clc

% test LMM for crossed designs, where the nested factor is null, but the
% coarser increase in level (x-axis)

%%% settings
nperm = 1000; % number of permutations in the ASCA model
%Slinspan = [0 exp(linspace(log(1e-4),log(0.0001),4))];
KK = 0:3;
niter = 1000; % for each setting, run niter models

% length(Slinspan)*niter;
DOE = fullfact([length(KK) niter]);
DOE(:,1) = KK(DOE(:,1));

% set seed
rng(321)

p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D = [D1(:,2) D1(:,1)]; % ID, condition, trt
lbD = {'ID','trt'};

% make unbalanced
% remove 50% of samples D2==1 and D1 in 1 and 2
% remove 50% of samples D2==2 and D1 in 4 and 5
% icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>1;
% icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>1;


% MODEL
M = [1 0; 0 1];
s = [0.2 0 0 .2];
%s(1) = Slinspan

plotit=0;
opt = ASCAcat('options');
opt.nperm = nperm;
opt.showtable = 0;

%%
%%% RECORD
pv = nan(size(DOE,1),2);
tic
for j=1:size(DOE,1)
    kk = DOE(j,1);
    % remove one at random
    ickeep = ones(size(D,1),1);
    if kk>0
        for ii=1:30
            rn = randperm(5);
            ickeep(D(:,1)==ii & ismember(D(:,2),rn(1:kk))) = 0;
        end
    end
    DD = D(ickeep==1,:);
    
    % simulate data
    [X DM Xm] = simulateData(DD,M,s,p,opt.contvar);
    
    % models 1 - Vanilla
    opt.rndfac = [];
    res1= ASCAcat(X,DD,M, opt);
    
    % models 2 - nested permutation
    opt.rndfac = 1;
    res2 = ASCAcat(X,DD,M, opt);
    
    %%% export results
    pv(j,:) = [res1.Effects{2}.p res2.Effects{2}.p];
    
    % save results
    if round(j / 100)==(j/100) || j < 3
        save('./simulations/BIGsimulation_crossedPerm.mat','DOE','pv')
    end
end

save('./simulations/BIGsimulation_crossedPerm.mat','DOE','pv')
toc

% %%
% clear PP
% for i=1:length(KK)
%     p = pv(DOE(:,1)==KK(i),:);
%     PP(i,:) = mean(p,1);
% end
% 
% plot(KK,PP);    shg

