clear; clc;
n = 30;
p = 40;
d1 = sort(repmat([1:2]',n/2,1));
x = randn(n,1)*3;
dep = [d1 x];
deplb = {'facA','covar'};
M = [1 0; 0 1];
s = [0.1 .1 0.1]
rndfac = []; 
deflate = 1;
contvar = 2; 
% generate data
Xc = sin(x*ones(1,p) + ones(n,1)*randn(1,p))*s(2);
Xd = mkdesignmatrix(d1)*randn(2,p)*s(1);
E = randn(n,p)*s(3);

X = Xc + Xd + E;
X = mncn(X); 

relD =  trace(Xd'*Xd)/ trace(X'*X);
relx =  trace(Xc'*Xc)/ trace(X'*X);
relsse =  trace(E'*E)/ trace(X'*X);

results = ASCAcat(X,dep,deplb,1000,M,rndfac,deflate,contvar);

results.ANOVAtab
[relD relx relsse]
%%
% make projection onto K
H = linspace(0.001,2,30)
for i = 1:length(H)
    [K h] = setCovariateKernel(x(randperm(n)),H(i));
    Xchat_crude = K*X; % crude estimate
    R(i,1) = corr(Xc(:),Xchat_crude(:));
    % marginal 1 (deflate by D1)
    D1 = mkdesignmatrix(d1);
    Xchat = K*(eye(n) - D1*pinv(D1))*X;
    R(i,2) = corr(Xc(:),Xchat(:));  
end

[K h] = setCovariateKernel(x);
% joint as sum of sequential deflated fits
Xdhat = D1*pinv(D1)*(X - Xchat);
Xhat = Xdhat + Xchat;
Ehat = X - Xhat;

subplot(2,3,1); plot(Xc(:),Xchat(:),'*'); abline(1,0); title('Marginal x'); shg
subplot(2,3,2); plot(Xd(:),Xdhat(:),'*'); abline(1,0); title('Marginal d'); shg
subplot(2,3,3); plot(X(:),Xhat(:),'*'); abline(1,0); title('all Hat'); shg
subplot(2,3,4); plot(E(:),Ehat(:),'*'); abline(1,0); title('resid'); shg
subplot(2,3,5); plot(Xc(:),Xchat_crude(:),'*'); title('Crude x'); abline(1,0); shg
[u ss v]= svds(Xchat,2);
subplot(2,3,6);
plot(x,u,'*'); shg
%%



%
% kernel smoothing regression
%n = length(x);
%h=median(abs(x-median(x)))/0.6745*(4/3/n)^0.2;
%kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
%rx=linspace(min(x),max(x),n*2);
%rx = sort(x);
%rx = x;
%K = [];
%for k=1:length(x)
%    z=kerf((x(k)-x)/h);
%    z = z / sum(z);
%    K(k,:) = z;
%f(k)=sum(z.*y);
%ff(k)=z'*y;
%end


y = sin(x) + randn(length(x),1)*0.1;
subplot(1,2,1)
plot(x,y,'*'); hold on
plot(x,K*y,'g*-');
%plot(rx,f+2,'*-');
%plot(rx,ff,'^-');
hold off; shg
subplot(1,2,2)
plot(y,K*y,'*'); shg
abline(1,0)