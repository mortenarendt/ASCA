function [Xs mx sx] = scale(X,mx,sx)

n = size(X,1); 
if nargin==1
    mx = nanmean(X); 
    sx = nanstd(X); 
end
Xs = X - (ones(n,1)*mx); 
Xs = Xs.*(ones(n,1)*(1./sx)); 


