function [D1t grps] = collapse(D1,unD1,TT)
%%
D1t = D1;

unTT = unique(TT);
grps = [];
for i=1:length(unTT)
    unD1sel = unD1(TT==unTT(i));
    if i==length(unTT)
        grps = [grps ['{' num2str(unD1sel') '}']];
    else
        grps = [grps ['{' num2str(unD1sel') '}, ']];
    end
    D1t(ismember(D1,unD1sel)) = unD1sel(1);
end