function m = makeScores(m)
% Function for calculating ASCA scores for each part of the design. 
% That is; Calculate loadings on the effect estimates, and calculate scores
% on the Effect estimates + E from those loadings. 

for i=1:length(m.Effects)
    XE = m.Effects{i}.Xd_crude + m.E; 
    totEV = sum(XE(:).^2);
    [u s v] = svds(m.Effects{i}.Xd_crude,rank(m.Effects{i}.Xd_crude));
    m.Effects{i}.loads{2} = v;
    
    Tm = m.Effects{i}.Xd_crude*v;
    Te = m.E*v;
    T = Tm + Te; 
    
    m.Effects{i}.loads{1} = T; 
    
    m.Effects{i}.varExp = [diag(T'*T) diag(Tm'*Tm) diag(Te'*Te)] / totEV;
    m.Effects{i}.varExp_description = {'colounm one is variation explained for Xm+E';...
        'colounm two is variation explained for Xm';...
        'colounm three is for E';'All relative to SSQ(Xm+E)'};
    
    m.Effects{i}.varExp2 = diag(Tm'*Tm) / sum(diag(Tm'*Tm));
    m.Effects{i}.varExp2_description = 'Variation explained for Xm relative to SSQ(Xm)';
end

