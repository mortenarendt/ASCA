function [X DM Xm] = simulateData(dep,M,s,p, contvar)
%s = [0.05 0 0 0.01];

nterms = size(M,1);
order = sum(M');
contterm = zeros(1,nterms);
for i=1:nterms
    nf = sum(M(i,:));
    iddep = find(M(i,:));
    contterm(i) = 0;
    if nf==1 % Main effects;
        % categorical variable
        if ismember(iddep,contvar)
            % continous variable
            DM{i} = setCovariateKernel(dep(:,iddep));
            contterm(i) = 1;
        else
            % categorical variable
            DM{i} = mkdesignmatrix(dep(:,iddep));
        end
    elseif nf==2 % two-way interaction effects.
        [DM{i},~,~,IDvec(:,i)] = mkdesignmatrix(dep(:,iddep(1)),dep(:,iddep(2)));
    elseif nf>2
        [DM{i},IDvec(:,i)] = mkdesignmatrixHigherOrder(dep(:,iddep));
    end
end


X = zeros(size(DM{1},1),p);
for i=1:size(M,1)
    iddep = find(M(i,:));
    if contterm(i)==1
        Xi = s(i)*getContinousHat(dep(:,iddep),p);
    else
        Bi = s(i)*randn(size(DM{i},2),p);
        Xi = DM{i}*Bi;
    end
    % center
    Xi = mycenter(Xi);
    Xm{i} = Xi;
    X = X+Xi; 
end
E = randn(size(X))*s(i+1);
E = mycenter(E); 
Xm{i+1} = E;
X = X + E;

function X = mycenter(X)
n = size(X,1); 
mn = mean(X,1);
X = X - ones(n,1)*mn; 
