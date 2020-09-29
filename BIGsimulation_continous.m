close all hidden; clear; clc
%%% settings
nperm = 0; % number of permutations in the ASCA model
Slinspan = [0 exp(linspace(log(1e-4),log(0.1),10))]; % loop over ev linear scalar
niter = 1000; % for each setting, run niter models

length(Slinspan)*niter;
DOE = fullfact([length(Slinspan) niter]);
DOE(:,1) = Slinspan(DOE(:,1));

% set seed
rng(123)

p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
covariate = randn(size(D1,1),1);
D = [D1(:,2) D2 D1(:,1) covariate]; % ID, condition, trt
lbD = {'ID','condition', 'trt','covar'};

% make unbalanced
% remove 50% of samples D2==1 and D1 in 1 and 2
% remove 50% of samples D2==2 and D1 in 4 and 5
icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>1;
icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>1;
icout = icout1 | icout2;
D = D(~icout,:);

% MODEL
M = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1;  0 1 1 0 ];
s = [0.1 0.1 0.1 0 0 .2];


plotit=0;
opt = ASCAcat('options');
opt.contvar = 4;
opt.contvaraslinear = 1;
opt.nperm = nperm;
opt.showtable = 0;


%%% RECORD
% EXPVAR for the three models including the oracle
% for selecting the correct model

EXPVAR = nan(size(DOE,1),4);
RELFIT  = nan(size(DOE,1),3);
pv = nan(size(DOE,1),3);
tic
for j=1:size(DOE,1)
    j
    % new covariate
    D(:,4) = randn(size(D,1),1);
    s(4) = DOE(j,1);
    % simulate data
    [X DM Xm] = simulateData(D,M,s,p,opt.contvar);
    
    % models 1
    opt.contvaraslinear = 1;
    opt.slopecorr = 0;
    reslin = ASCAcat(X,D,M, opt);
    
    % models 2
    opt.contvaraslinear = 0;
    opt.slopecorr = 0;
    ressm = ASCAcat(X,D,M, opt);
    
    % models 3
    opt.slopecorr = 1;
    ressmslc = ASCAcat(X,D,M, opt);
    
    %%% export results
    EXPVAR(j,:) = [trace(Xm{4}'*Xm{4}) reslin.Effects{4}.ExpVar_D_marginal ressm.Effects{4}.ExpVar_D_marginal ressmslc.Effects{4}.ExpVar_D_marginal];
    
    df = Xm{4} - reslin.Effects{4}.Xd_marginal;
    
    r1 = corr(Xm{4}(:),reslin.Effects{4}.Xd_marginal(:));
    RELFIT(j,1) = trace(df'*df) / trace(Xm{4}'*Xm{4});
    df = Xm{4} - ressm.Effects{4}.Xd_marginal;
    r2 = corr(Xm{4}(:),ressm.Effects{4}.Xd_marginal(:));
    RELFIT(j,2) = trace(df'*df) / trace(Xm{4}'*Xm{4});
    df = Xm{4} - ressmslc.Effects{4}.Xd_marginal;
    r3= corr(Xm{4}(:),ressmslc.Effects{4}.Xd_marginal(:));
    RELFIT(j,3) = trace(df'*df) / trace(Xm{4}'*Xm{4});
    CORR(j,:) = [r1 r2 r3 ]; 
    pv(j,:) = [reslin.Effects{4}.p ressm.Effects{4}.p ressmslc.Effects{4}.p]; 
    S2N(j) = trace(Xm{4}'*Xm{4}) / trace(Xm{end}'*Xm{end});
    XX = Xm{1} + Xm{2} + Xm{3} + Xm{5};
    background2N(j) = trace(XX'*XX) / trace(Xm{end}'*Xm{end});
    % save results 
end
toc
%%
save('./simulations/BIGsimulation_continous3.mat','EXPVAR','DOE', 'S2N','RELFIT', 'CORR')
