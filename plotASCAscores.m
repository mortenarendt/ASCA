function plotASCAscores(m,eff,cmp)

% plots ASCA scoreplot 
% input: 
%   m = ASCA model (output from ASCAcat())
%   eff = which model parameter to plot 
%   cmp = which components - default 1 vs 2

if nargin==2
    cmp = [1 2];
end

u = m.Effects{eff}.loads{1};
EV = m.Effects{eff}.varExp; 
cl = m.dependent_var(:,eff);
tit = m.dependent_var_lb{eff}; 
EVeff = str2num(cell2mat(m.ANOVAtab(ismember(m.ANOVAtab(:,1),tit),3)))*100;

if cmp(2) > size(u,2)
    % add a vector from the residual matrix
    [uu ss vv] = svds(m.E,1);
    EV = [EV; [0 0 0]];
    u = [u uu]; 
end

scatter(u(:,cmp(1)),u(:,cmp(2)),100,cl,'filled'); 

title([tit ' (' num2str(EVeff) '% - of Total Variance)']); 
xlabel(['Component ' num2str(cmp(1)) ' - ' num2str(EV(cmp(1),1)*100,3) ...
    '% (model: ' num2str(EV(cmp(1),2)*100,3)...
    '%  - resid: ' num2str(EV(cmp(1),3)*100,3) '%)'])

ylabel(['Component ' num2str(cmp(2)) ' - ' num2str(EV(cmp(2),1)*100,3) ...
    '% (model: ' num2str(EV(cmp(2),2)*100,3)...
    '%  - resid: ' num2str(EV(cmp(2),3)*100,3) '%)'])

