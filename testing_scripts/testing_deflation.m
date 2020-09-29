clear; clc;
nrep = 6;
nA = 5;
nB = 2;
p = 40;
D = repmat(fullfact([nA,nB]),nrep,1);

N = size(D,1);
q = 0;
% make unbalanced
icout1 = D(:,2)==1 & ismember(D(:,1),[1 2]) & rand(N,1)<q;
icout2 = D(:,2)==2 & ismember(D(:,1),[4 5]) & rand(N,1)<q;
icout = icout1 | icout2;
D = D(~icout,:);
n = size(D,1); 
x = randn(n,1)*3;
dep = D;
deplb = {'FacA','Facb'};
M = [1 0; 0 1; 1 1];
rndfac = [];
deflate = 1;
contvar = [];

s = [0.1 .1 0.2 0.0001]
% generate data
X = zeros(size(D,1),p);
Djoint  = []; 
for i=1:2
    Bi = s(i)*randn(length(unique(D(:,i))),p);
    Xi = mkdesignmatrix(D(:,i))*Bi;
    X = X+Xi;
    
    d = mkdesignmatrix(D(:,i)); 
    
    PX{i} = d*pinv(d);
    Djoint = [Djoint d]; 
end

Xc = sin(x*ones(1,p) + ones(n,1)*randn(1,p))*s(i+1);
X = X + Xc + randn(size(X))*s(i+2);
X = mncn(X);
%results = ASCAcat(X,dep,deplb,1000,M,rndfac,deflate,contvar)


for i=3:4
    dd = ceil(rand(n,1)*5);
    d = mkdesignmatrix(dd); 
    dep = [dep dd];
    PX{i} = d*pinv(d);
    Djoint = [Djoint d]; 
end
i = 5; 
K = setCovariateKernel(x);
dep = [dep x];
PX{i} = K;

%deplb = [deplb {'FacN1'} {'FacN2'} {'cov'}];
%M = [eye(5); 1 1 0 0 0]; 
dep = dep(:,1:2); 
M = [1 0; 0 1]
results = ASCAcat(X,dep,deplb,1000,M,rndfac,deflate,5);
%%
results.ANOVAtab
%%
clc
[P1 P2] = blockwisehat(PX([1 2 3 4 5]));

%%% check it
E1 = (eye(n) - P1 - P2)*X; 
%E2 = (eye(n) - P1)*X; 
%plot(E1(:),E2(:),'*'); shg
%trace(E1'*E1)
%trace(E2'*E2)

sqrt(trace(E1'*E1) / ((n-1)*(p)))
s(4)


subplot(1,2,1); plot(x,PX{5}*X(:,2),'*'); shg
subplot(1,2,2); plot(x,P2*X(:,2),'*'); shg

%%
n = size(D,1);
D1 = mkdesignmatrix(D(:,1));
D2 = mkdesignmatrix(D(:,2));

W1 = D1*pinv(D1);
W2 = D2*pinv(D2);
W1 = setCovariateKernel(x);

DefM1 = (eye(n) - W1);
DefM2 = (eye(n) - W2);
%W2o1 = (DefM1*D2)*pinv(DefM1*D2);
%W2o1 = DefM1*D2*pinv(DefM1*D2);


%DefM2o1 = (eye(n) - W2o1);
%E1 = DefM1*DefM2*X;
%E1 = DefM2*DefM1*DefM2*X;
%E1 = (eye(n) - DefM1*W2*DefM1')*DefM1*X;

%PX = PA + P(MAB)
%PA = W1
%P(MAB) = (DefM1*W2)*pinv(DefM1*W2)
%PX = W1 + (DefM1*W2)*pinv(DefM1*W2);

PX2 = W2;
PX1 = (DefM2*W1)*pinv(DefM2*W1);
%PX = PX2 + PX1;
E1 = (eye(n) - PX2 - PX1)*X;
SSEm = trace(E1'*E1);
%E1 = (eye(n) - W2o1)*DefM1*X;

% scrample PX1
for i=1:100
    idd = randperm(n);
    PX1p = (DefM2*W1(idd,idd))*pinv(DefM2*W1(idd,idd));
    E1p = (eye(n) - PX2 - PX1p)*X;
    SSEp(i) = trace(E1p'*E1p);
end
%D2m = DefM1*D2;
%W2m = D2m *pinv(D2m);
%E1 = (eye(n) - W2m)*DefM1*X;

%E2 = (eye(n) - [D1 D2]*pinv([D1 D2]))*X;
%trace(E2'*E2)

%hist(SSEp); vline(SSEm); 
pv = (sum(SSEp<SSEm)+ 1) / length(SSEp)
%shg






