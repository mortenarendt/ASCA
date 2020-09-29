clear; clc;
n = 30;
p = 40;
niter = 100;
deplb = {'facA','covar'};
M = [1 0; 0 1];
s = [0.1 0 0.1]
rndfac = [];
deflate = 1;
contvar = 2;
R = 0:0.1:1;
d1 = sort(repmat([1:2]',n/2,1));
c = 0; 
for j=8 %1:length(R)
    r = R(j);
    r
    for ii=1 %:niter
        c = c+1;
        x = r*(d1-1.5) + (1-r)*randn(n,1);
        dep = [d1 x];
        
        % generate data
        Xc = sin(x*ones(1,p) + ones(n,1)*randn(1,p))*s(2);
        Xd = mkdesignmatrix(d1)*randn(2,p)*s(1);
        E = randn(n,p)*s(3);
        
        X = Xc + Xd + E;
        X = mncn(X);
       
        
        results = ASCAcat(X,dep,deplb,1000,M,rndfac,deflate,contvar);
        pv(c,:) = [r ii results.Effects{2}.p_F results.Effects{2}.p_F_est];
    end
end
%%       
RES = []; 
for j=1:length(R)
    RES = [RES; mean(pv(pv(:,1)==R(j),:),1)]; 
end

plot(RES(:,1),RES(:,3), '*-'); shg

