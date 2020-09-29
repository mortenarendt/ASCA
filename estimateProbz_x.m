function Pz_x = estimateProbz_x(X,Z)

% establish conditional distribution distribution Q(X | Z)
% fit line and uncertainty
n = length(X);
b = pinv([ones(n,1) X])*Z;
zhat = [ones(n,1) X]*b;
e = Z - zhat;
se = sqrt(e'*e / (n-1));


%% for every X predict likelihood of all Z
%Zdiff = Z*ones(1,n) - ones(n,1)*zhat';
Zdiff = zhat*ones(1,n) - ones(n,1)*zhat';
%Pz_x = 1 - normpdf(abs(Zdiff) / (se));
Pz_x = 1 - normcdf(abs(Zdiff) / (se*sqrt(2/n)));
Pz_x  = normaliz(Pz_x);

Zdiff2 = (zhat/se)*ones(1,n) - ones(n,1)*(zhat'/se);
%Pz_x2 = 1 - normcdf(abs(Zdiff2) / (1*sqrt(2/n)));
Pz_x2 = 1 - chi2cdf(Zdiff2.^2 ,2);
Pz_x2  = normaliz(Pz_x2);
plot(Pz_x(:),Pz_x2(:),'.'); shg
%%
x = 0:0.01:30;
pv = chi2pdf(x,2); 

z = randn(10000,1) - randn(10000,1); 
[h p] = hist(z.^2,50); 
pp = p(2) - p(1);
xx = 0.01;
bar(p,h / sum(h*pp)); hold on; 
plot(x,pv / sum(pv*xx),'r-', 'linewidth',2); hold off; shg




function Xn = normaliz(X)
n = size(X,2);
Xn  = X./(sum(X,2)*ones(1,n));
