function id = nestedrandperm(nestedfac,classD)
% id = nestedrandperm(nestedfac, classD)
% return a randomsation of the classD labels but such that these are
% contained nested within nestedfac (which should be a finner version of
% classD)

%%
% make new cross table btw factor of interest and the nested version
[un, ID1,ID] = unique(nestedfac);

%for jj=1:100
unp = un(randperm(length(un)));
cl2 = classD(ID1);
uncl2 = unique(classD);
% find wihtin each group which nested factors should go in
id = zeros(1,length(classD));

for j=1:length(unp)
    candiateIDs = find(ID==unp(j));
    nspots = sum(ID==un(j));
    
    % if needed repeat some of the observations in the permutation
    cid = [];
    for i = 1:ceil(    nspots / length(candiateIDs ))
        cid = [cid candiateIDs(randperm(length(candiateIDs)))];
    end
    % reduce to the number of spots in the id vector
    cid = cid(1:nspots);
    
    id(ID==un(j)) = cid;
end

%%
%check
%
classDp = classD(id);

%for j=1:length(un)
%    CL(jj,j) =unique(classDp(nestedfac==un(j)));
%    length(unique(classDp(nestedfac==un(j))))
%end

%end
%%
%sum(CL==1)
%[nestedfac classDp]