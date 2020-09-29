close all hidden; clear; clc
p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
cv = rand(size(D1,1),1);
D = [D2 D1(:,1) cv];
lbD = {'condition', 'trt','cov'};

s = [0 0 0.1 0.01];
M = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
rndfac = [];

X = zeros(size(D,1),p);
for i=1:2
    Bi = s(i)*randn(length(unique(D(:,i))),p);
    Xi = mkdesignmatrix(D(:,i))*Bi;
    X = X+Xi;
end
i = i+1;
Bi = s(i)*randn(1,p);
Xi = D(:,i)*Bi;
X = X+Xi;
X = X + randn(size(X))*s(4);
%%
%X = mncn(X);
res = ASCAcat(X,D,lbD,100,M,rndfac,1,[3]);