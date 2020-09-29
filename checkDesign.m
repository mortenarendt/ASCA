function tab = checkDesign(D,Dlb)

% investigate the structure of the design by investigation of the number of
% replicates within each design-cell. 
[n p] = size(D);

T = [];
for i=1:p
    d = D(:,i);
    if i<p
        d  = [ repmat([char(Dlb(i)) '='],n,1) num2str(d) repmat(', ',n,1)];
    else
        d  = [ repmat([char(Dlb(i)) '='],n,1) num2str(d)];
    end
    T = [T d];
    
end
tab.T = cellstr(T);
[tab.undesign idu,ida] = unique(tab.T);
tab.h = hist(ida,unique(ida));
barh(1:length(tab.h),  tab.h);
text(tab.h,1:length(tab.h), tab.undesign);
xlim([0 max(tab.h) + 3])
