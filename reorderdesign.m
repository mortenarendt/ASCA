function ordr = reorderdesign(PX)

% figure out the order of the first nterms-1 
nterms = length(PX); 
prms = perms(1:(nterms));
for i=1:size(prms ,1)
    px = zeros(size(PX{1})); 
    for ii=1:size(prms,2)
        px = px + PX{prms(i,ii)};
        rnk(i,ii) = rank(px);
    end
end


%%% find an increasing combination
Drnk = rnk(:,2:end) - rnk(:,1:end-1);
ic1 = sum(Drnk>0,2)==(size(rnk,2)-1);
rnk = rnk(ic1,:);
prms = prms(ic1,:);
%%% choose the one closest to the input (1:nterms-1)
[~,id] = max(prms*[1:nterms]');
ordr = prms(id,:); 