close all hidden; clear; clc
nperm = 100;
nsim = 100;
p = 30; % number of variables
n = 10;
s = [0 0.1 0.01]; % effect size of factors
M = [1 0; 0 1];
nA = 2;
nB = 5;
lbD = {'A', 'B'};

pv1 = [];
pv2 = [];
qq = [];

Q = 0:0.1:1;
for ii=1:length(Q)
    q = Q(ii)
    % calculate nrep.
    % 6*nrep + 4*nrep*(1-q) = nrep*(6 + 4(1-q)) = 10*n: => nrep = 10*n / (6 +
    % 4*(1-q))
    nrep = round(10*n / (6 + 4*(1-q)));
    
    D = repmat(fullfact([nB nA]),nrep,1);%
    N = size(D,1);
    % make unbalanced
    icout1 = D(:,2)==1 & ismember(D(:,1),[1 2]) & rand(N,1)<q;
    icout2 = D(:,2)==2 & ismember(D(:,1),[4 5]) & rand(N,1)<q;
    icout = icout1 | icout2;
    D = D(~icout,:);
    
    %%% for this design, run 100 simulations
    for j = 1:nsim
        % generate data
        X = zeros(size(D,1),p);
        for i=1:2
            Bi = s(i)*randn(length(unique(D(:,i))),p);
            Xi = mkdesignmatrix(D(:,i))*Bi;
            X = X+Xi;
        end
        X = X + randn(size(X))*s(i+1);
        X = mncn(X);
        res_def = ASCAcat(X,D,lbD,nperm,M,[],1,[]);
        res_nodef = ASCAcat(X,D,lbD,nperm,M,[],0,[]);
        
        % take out pv for design A.
        pv1 = [pv1 res_def.Effects{1}.p];
        pv2 = [pv2 res_nodef.Effects{1}.p];
        qq = [qq q];
    end
end
%%
save('simulate2fac_balanced_vs_unbalanced_v2.mat','pv1','pv2','qq')
%%
close all;
ntest = size(pv,2);
for i=1:ntest
    subplot(2,ntest,i);
    hist(pv(:,i),30);xlim([0 1])
    subplot(2,ntest,i+ntest);
    hist(pv_est(:,i),30);xlim([0 1])
end
shg;



