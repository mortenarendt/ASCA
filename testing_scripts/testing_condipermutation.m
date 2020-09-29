clear; 
clc;
%close all;
niter = 100;
nrep = 30;
RES = [];
%Sample data for different size
N = [20 50 100 200 500];
N = [100 200 500];
for j=1:length(N)
    n = N(j)
    for jj = 1:nrep
        X = sort(randn(n,1));
        Z = randn(1)*X + randn(1)*randn(n,1);
        Pz_x = estimateProbz_x(X,Z);
        
        cr = [];
        for i = 1:niter
            id = permutesampler(Pz_x);
            cr(i) = corr(X(id),Z);
        end
        
        crobs = corr(X,Z);
        RES = [RES; [n crobs mean(cr - crobs) median(cr - crobs)]];
    end
end


save testing_conditionalpermutation.mat RES niter nrep
%%
cr = [];
for i = 1:niter
    id = condirandomperm(X,Z);
    
end

hist(cr,30); xlim([-1,1])
vline(corr(X,Z),'r');
vline(mean(cr),'b'); shg

%%

% establish conditional distribution distribution Q(X | Z)
% fit line and uncertainty
b = pinv([ones(n,1) X])*Z;
e = Z - [ones(n,1) X]*b;
se = sqrt(e'*e / (n-1));
% for every X predict likelihood of all Z
zhat = [ones(n,1) X]*b;
Pz_x = normpdf(abs((Z*ones(1,n) - ones(n,1)*zhat') / se));
Punif = rand(n,n);

smpl = Punif.*Pz_x;

id = [];
rorder = randperm(n);
for i=rorder
    % for each row take the most likeli observation and remove it from the
    % possible set
    [~,id(i)] = max(smpl(i,:));
    smpl(:,id(i)) = 0;
    
end
length(unique(id))
id
contourf(Punif.*Pz_x); shg

%%
clear; clc;
nrep = 3;
nA = 5;
nB = 2;
D = repmat(fullfact([nA,nB]),nrep,1);
N = size(D,1);
q = 0.5;
% make unbalanced
icout1 = D(:,2)==1 & ismember(D(:,1),[1 2]) & rand(N,1)<q;
icout2 = D(:,2)==2 & ismember(D(:,1),[4 5]) & rand(N,1)<q;
icout = icout1 | icout2;
%D = D(~icout,:);

% establish conditional distribution
%Pij = P(D1=d1i | D2==d2j)
%%
p = [0.25 0.5 0.25];
for i=1:1000
    %[~,id(i)] = min(abs(cumsum(p) - rand(1)));
    id(i) = sum((cumsum(p) - rand(1))<0)+1;
end
hist(id); shg




