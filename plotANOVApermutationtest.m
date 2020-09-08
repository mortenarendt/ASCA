function plotANOVApermutationtest(res)

np = length(res.Fperm); 
% hist(ss3fp); vline(sse3f); shg;
[h, pp] = hist(res.Fperm,round(sqrt(np*2)));
bar(pp,h);
hh = vline(res.Fmodel); hold on; shg;
set(hh,'linewidth',2)
xlm = get(gca,'xlim');
xax = linspace(xlm(1),xlm(2)+(xlm(2)-xlm(1))*0.05,1000);
plot(xax,fpdf(xax,res.dftop_est,res.dfbottom_est)*res.nperm*(pp(2)-pp(1)),'r-','linewidth',2)
plot(xax,fpdf(xax,res.dftop,res.dfbottom)*res.nperm*(pp(2)-pp(1)),'c-','linewidth',2)
hold off;
axis tight;
xlm = get(gca,'xlim');
ylm = get(gca,'ylim');
xpt = xlm(1) + 0.05*(xlm(2)-xlm(1));
xpt2 = xlm(1) + 0.25*(xlm(2)-xlm(1));
ypt = ylm(2)*0.9;
text(xpt,ypt,['p freq = ' num2str(res.p,3)]);
text(xpt2,ypt*0.8,['p est = ' num2str(res.p_F_est,3)]);
xlabel('F');
try
title([res.label [' DF = ' num2str(res.rank)]],'interpret','none')
end