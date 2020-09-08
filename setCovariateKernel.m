function [K h] = setCovariateKernel(x,h)
% K = setCovariateKernel(x,h)
% x ~ a n,1 contious vector
% h ~ bandwith (can be omitted)


n = length(x);
x = x(:);

if nargin==1
    h=median(abs(x-median(x)))/0.6745*(4/3/n)^0.2;
end


kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
K = [];
for k=1:length(x)
    z=kerf((x(k)-x)/h);
    z = z / sum(z); % normalize to unit norm
    K(k,:) = z;
end



