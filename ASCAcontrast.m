function results = ASCAcontrast(res1, which_factor)

% Contrast estimation based on an ASCA model output. 
% Using the marginal effect matrix associated with the factor
% (which_factor), First a hieracical cluster is produced. At each hight of
% this dendrogram, the levels below are lumped together, and a marginal
% test is run towards the full model. 
% I/O: 
% Input: res1 = ASCAcat() output, which_factor = the factor of
% interest
% Output: results: dendrogram, and test at each level. 
% plot the results by: plotASCAcontrast(results) 

results.originalmodel = res1; 
results.which_factor = which_factor; 
M = res1.model;
D = res1.dependent_var; 
X = res1.data; 

nperm = zeros(1,size(M,1)+1);
nperm(which_factor) = res1.options.nperm(which_factor);
opt2 = res1.options;
opt2.nperm = nperm;

% contrast estimation for factor 1
% 1) Take effect matrix and get dendrogram for comparison

Xd = res1.Effects{which_factor}.Xd_marginal;
D1 = res1.dependent_var(:,which_factor);

[unD1, idd] = unique(D1);

results.Xd = Xd(idd,:);
results.groups = unD1; 

results.distY = pdist(results.Xd);


results.linkageZ = linkage(results.distY);

nf = length(idd);
T = cluster(results.linkageZ,'maxclust',(nf-1):-1:2);

RES = []; 
for i=1:size(T,2)
    [D1t grps{i}] = collapse(D1,unD1,T(:,i));
    
    % 2) run marginal models with a collapsed design on the individual parts.
    % example collapse d [1 2 3] to one group and run model with marginal
    
    DD = [D D1t];
    MM = [M zeros(size(M,1),1); zeros(1,size(M,2)) 1];
    
    
    res1t= ASCAcat(X,DD,MM, opt2);
    
    % Record p and F
    RES(i,:) = [length(unique(D1t))  res1t.Effects{which_factor}.Fmodel res1t.Effects{which_factor}.p]; 
end

results.hierachical_contrasts = RES;
results.hierachical_contrasts_groups = grps;
