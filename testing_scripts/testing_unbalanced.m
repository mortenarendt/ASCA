close all hidden; clear; clc
p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
D = [D1(:,2) D2 D1(:,1) randn(size(D1,1),1)]; % ID, condition, trt
lbD = {'ID','FacA', 'FacB','cov'};
%
% make unbalanced
% remove 50% of samples D2==1 and D1 in 1 and 2
% remove 50% of samples D2==2 and D1 in 4 and 5
icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>0.5;
icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>0.5;
icout = icout1 | icout2;
sum(icout)
D = D(~icout,2:end);
lbD = lbD(2:end);
n = size(D,1); 

M = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
s = [1 1 1 1];

rndfac = [];
X = zeros(size(D,1),p);
% main effects ONLY
for i=1:2
    Bi = s(i)*randn(length(unique(D(:,i))),p);
    Xi = mkdesignmatrix(D(:,i))*Bi;
    X = X+Xi;
end

Xc = sin(D(:,3)*ones(1,p) + ones(n,1)*randn(1,p))*s(3);
X = X + Xc + randn(size(X))*s(4);

X = mncn(X);
%
%results = ASCAcat(X,dep,deplb,nperm,M,rndfac,deflate,contvar)
res = ASCAcat(X,D,lbD,100,M,rndfac,1,3);