function [D, d1label, d2label, class] = mkdesignmatrix(d1,d2,nan,verbose)

if nargin==1 %d1
   nan = 1;  
    [D d1label] = mkmatrix(d1,nan);
    return
    
elseif nargin==2
    if numel(d2)==1 %d1 and nan;
        nan=d2;
        [D d1label] = mkmatrix(d1,nan);
        return
    else % d1 and d2
        nan = 1;
    end
end

if nargin<4
    verbose = 0; 
end
if verbose==1
disp('Two input design matrixes - interaction matrix is produced')
end

D1 = mkmatrix(d1,nan);
D2 = mkmatrix(d2,nan);

D1ex = repmat(D1,1,size(D2,2));
d1label = repmat(1:size(D1,2),1,size(D2,2));
id = sort(repmat(1:size(D2,2),1,size(D1,2)));
D2ex= D2(:,id);
d2label = id;
sumD = D1ex + D2ex;
D = sumD==2;
D = D+0;

[r c] = find(D);
[sr id] = sort(r);
class = zeros(size(D,1),1);
class(sr) = c(id);



function [D du] = mkmatrix(d,nan)
if isnumeric(d)
    du = unique(d);
    idnan = find(isnan(du));
    if ~isempty(idnan);
        if nan==1;
            du = du(1:idnan(1));
        elseif nan==0;
            du = du(isnan(du)==0);
        end
    end
    
    D = zeros(length(d),length(du));
    for i=1:length(du);
        if isnan(du(i))
            D(isnan(d),i) = 1;
            D(isnan(d),1:i-1) = NaN;
        else
            D(d==du(i),i) = 1;
        end
    end
    
elseif iscell(d)
    du = unique(d);
    D = zeros(length(d),length(du));
    for i=1:length(du);
        if ismember(du(i),cellstr(''));
            idnan = ismember(d,'');
            D(idnan,i) = 1;
            D(idnan,1:i-1) = NaN;
        else
            D(ismember(d,du(i)),i) = 1;
        end
    end
end

