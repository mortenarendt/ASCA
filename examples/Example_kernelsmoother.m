close all hidden; clear; clc
rng(15)
n = 30; 
u = sort(rand(n,1) - 0.5); 
u2 = linspace(min(u),max(u),1000)';
X(:,1) = sin(5*[u; u2]);
ord = 3; 
B = [ones(length([u; u2]),1)];
for i=1:ord
    B(:,i+1) = [u; u2].^i;
end
X(:,2) = B*randn(ord+1,1);
X(:,3) = 1./ ([u; u2] - min(u) + 0.1);

% scale to the same min / max
X = X - ones(size(X,1),1)*min(X); 
X = X*diag(1./max(X)); 
E = randn(n,3)*0.1;

% kernel fit
K = setCovariateKernel(u,0.07);
Xh = K*(X(1:n,:)+E);

X2 = X((n+1):end,:);
X = X(1:n,:);
% plot it
plot(u2,X2,'-'); hold on; 
plot(u,X+E,'*'); 
plot(u,Xh,'^'); hold off; 

save kernelsmoother_example.mat




