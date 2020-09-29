function id = condirandomperm(X,Z)

% establish conditional distribution distribution Q(X | Z)
% fit line and uncertainty
n = length(X);
b = pinv([ones(n,1) X])*Z;
zhat = [ones(n,1) X]*b;
e = Z - zhat;
r = 1 - (e'*e) / (Z'*Z); 
r = sqrt(r); 

se = sqrt(e'*e / (n-1));

% for every X predict likelihood of all Z
%Zdiff = Z*ones(1,n) - ones(n,1)*zhat';
Zdiff = zhat*ones(1,n) - ones(n,1)*zhat';
%Pz_x = 1 - normpdf(abs(Zdiff) / (se));
Pz_x = 1 - normcdf(abs(Zdiff) / (se*sqrt(2/n)));

Pz_x  = normaliz(Pz_x);
Punif = normaliz(rand(n)); 

smpl = r*Pz_x + (1-r)*Punif;
%smpl = Punif.*Pz_x;
%smpl = Punif + Pz_x;
%smpl = Pz_x;
%smpl = Punif; 

id = zeros(n,1);
%id2 = zeros(n,1);
rorder = randperm(n);
%rorder = 1:n;
idd = 1:n; 
for i=rorder
    % for each row take the most likeli observation and remove it from the
    % possible set
    %[~,id(i)] = max(smpl(i,:));
    %smpl(:,id(i)) = 0;  
    
    % sample for each row
    id2 = sum((cumsum(Pz_x(i,:)) - rand(1))<0)+1;
    id(i) = idd(id2); 
    
    % remove the coloumn
    Pz_x(:,id2) = [];
    idd(id2) = []; 
    Pz_x = normaliz(Pz_x);
    
    
    
    
end
%length(unique(id))
%id
%contourf(Punif.*Pz_x); shg

function Xn = normaliz(X)
n = size(X,2);
Xn  = X./(sum(X,2)*ones(1,n));

