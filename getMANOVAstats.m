function m= getMANOVAstats(H,E)

% get several stats 
%H = E2f - E3f; H = H'*H; 
%E = EtE;
%%% Phillais trace
philtrace = trace(H/(H + E));
%%% HotellingsLawley trace
hottrace = trace(H/E);
%%% Wilks Lambda
%%
ev = svds(H/E,rank(H));
%wilksL = prod(1./ev+1);
wilksL = det(E) / det(H + E);
%%
%%% Roys largest root
roysmaxev = max(ev);
m = [philtrace hottrace wilksL roysmaxev];
