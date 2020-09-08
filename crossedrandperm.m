function id = crossedrandperm(nestedfac,classD)
% id = crossedrandperm(nestedfac, classD)
% return a randomsation of the classD labels within groups indicated in nestedfac

%%
% make new cross table btw factor of interest and the nested version
[un, ~ ,ID] = unique(nestedfac);

% find wihtin each group which nested factors should go in
id = zeros(1,length(classD));
for j=1:length(un) 
    idd = find(ID==un(j)); 
    id(idd) = idd(randperm(length(idd))); 
end