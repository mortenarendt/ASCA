close all hidden; clear; clc
% asca setting
opt = ASCAcat('options');
opt.showtable = 0;
opt.rndfac = 1; % the first factor is random
opt.nperm = 50;

% set data and models variance contribution for each factor
p = 20;
M = [1 0 0; 0 1 0; 0 0 1; 0 1 1];
s = [0.05 0 0 0 0.01];

% set design parameters
k = 30;
rep = 5;
THR = 0.1:0.1:1;
%%
for jj = length(THR)
    thr = THR(jj);
    jj
    D = getDesign(k,rep,thr);
    for j = 1:100
        X = simulateData(D,M,s,p,[]);
        res = ASCAcat(X,D,M, opt);
        RESULTS{j} = res;
        
        for i = 1:length(res.Effects)
            pv(j,i) = res.Effects{i}.p;
            pv_est(j,i) = res.Effects{i}.p_F_est;
            
        end
        
    end
    MN(jj,1,:) = median(pv);
    MN(jj,2,:) = mean(pv);
end


%%
plot(THR,squeeze(MN(:,2,:))); shg
%%
close all hidden;
close all;
ntest = size(pv,2);
for i=1:ntest
    subplot(2,ntest,i);
    hist(pv(:,i),30);xlim([0 1])
    subplot(2,ntest,i+ntest);
    hist(pv_est(:,i),30);xlim([0 1])
end
shg;



