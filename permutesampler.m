function id = permutesampler(Pz_x)

n = size(Pz_x,1);
id = zeros(n,1);
%id2 = zeros(n,1);
rorder = randperm(n);
%rorder = 1:n;
idd = 1:n; 
for i=rorder(1:end-1)
    % sample for each row
    id2 = sum((cumsum(Pz_x(i,:)) - rand(1))<0)+1;
    id(i) = idd(id2); 
    
    % remove the coloumn
    Pz_x(:,id2) = [];
    idd(id2) = []; 
    %Pz_x = normaliz(Pz_x);
    Pz_x = TSnorm(Pz_x);
    
end

id(rorder(end)) = idd; 

function Xn = normaliz(X)
n = size(X,2);
Xn  = X./(sum(X,2)*ones(1,n));