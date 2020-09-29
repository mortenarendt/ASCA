close all hidden; clear; clc

nperm = 1000;
niter = 100;

p = 30;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
lbD = {'ID','FacA', 'FacB'}; %,'cov'};
lbD = lbD(2:end)


s = [1 0 1 1];

rndfac = [];
M = [1 0; 0 1];

Q = 0:0.05:0.99
%
RES = [];lbs = {'Q','N','iter','propmatch','pv'};
for jjj=1%:niter
    jjj
    for jj=20 %:length(Q)
        
        q = Q(jj);
        D = [D1(:,2) D2 D1(:,1)]; % randn(size(D1,1),1)]; % ID, condition, trt
        % make unbalanced
        icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)<q;
        icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)<q;
        icout = icout1 | icout2;
        %sum(icout)
        D = D(~icout,2:end);
        n = size(D,1);
        
        X = zeros(size(D,1),p);
        % main effects ONLY
        for i=1:2
            Bi = s(i)*randn(length(unique(D(:,i))),p);
            Xi = mkdesignmatrix(D(:,i))*Bi;
            X = X+Xi;
        end
        
        %Xc = sin(D(:,3)*ones(1,p) + ones(n,1)*randn(1,p))*s(3);
        X = X + randn(size(X))*s(4);
        X = mncn(X);
        
        res1 = ASCAcat(X,D,lbD,nperm,M,rndfac,1,3,1);
        RES = [RES; [q n jjj 0 res1.Effects{1}.p]];
        %res1.Effects{1}.p
        
        %%% manual way
        d1 = mkdesignmatrix(D(:,2));
        d2 = mkdesignmatrix(D(:,1));
        for i=1:nperm
            idd = randperm(n); 
            E = (eye(n) - [d1(idd,:) d2]*pinv([d1(idd,:) d2]))*X; 
            ssep(i) = trace(E'*E);
        end
        
        %res2 = ASCAcat(X,D,lbD,nperm,M,rndfac,1,3,1);
        %res2.Effects{1}.p
        %RES = [RES; [q n jjj 1 res2.Effects{1}.p]];
    end
    %save testing_propmatching_woOrth.mat RES lbs
end
%
SSEp = [res1.Effects{2}.permutationSSE' ssep']; 
%%
clc
subplot(2,1,1); hist(SSEp(:,1),100); vline(res1.Effects{2}.modelSSE); shg    
subplot(2,1,2); hist(SSEp(:,2),100); vline(res1.Effects{2}.modelSSE); shg    
mean(SSEp,1)
sum(SSEp<res1.Effects{2}.modelSSE,1)/nperm

