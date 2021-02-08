function res = ANOVApermutationtest(x,PX,D,nperm,contvar,prop_matching,classD,nestedfac,propMcov )

% x - Data
% PX is a list of projection matrices, one for each term in the entire anova model
% The last element of PX is the factor of interest.
% D The design matrix similar to PX in organization
% contvar [0 or 1] indicates whether the factor of interest is computed
% from a continous predictor (1) - and hence a smoothing kernel, or a
% factor (0)
% prop_mathcing [0 or 1] indicates whether to use propensity score
% mathching in the scrambling procedure
% classD (vector) is the factor of interest
% nestedfac (vector) is the random term towards which the factor should be
% tested.

alpha_cut = 0.5;
plotit = 0;
[n p] = size(x);

% check if nested design (based on classD and nestedfactor)
if nargin<7
    % Normal ANOVA type model, NO mixed effects
    isnested = 0;
    iscrossed = 0;
else
    d1 = mkdesignmatrix(classD);
    d2 = mkdesignmatrix(nestedfac);
    isnested = (rank(d2)==rank([d1 d2]))+0;
    
    if isempty(d2) 
        isnested = 0; 
    end
    % is it crossed? or only nested
    if isnested==0
        iscrossed = 1;
    else
        iscrossed = 0;
    end
    if isempty(d1)
        isnested=0; iscrossed=0;
    end
end


[P1 P2] = blockwisehat2(D); % returns P2 - the factor of interest (last element of PX) and P1 the rest.
%E = (eye(n) - P1 - P2)*x;
%[P1n P2n] = blockwisehat2(D([2 1 3]));
%E2 = (eye(n) - P1n - P2n)*x;
%plot(E(:),E2(:),'*'); shg
%%
%P1 = PX{1}; P2 = PX{2};
defM1 = eye(n) - P1;
%df3f = rank(P1+P2)
%df2f = rank(D2);
%df2f = rank(P1);

df3f = trace(P1 + P2);
df2f = trace(P1);

% crude estimate (hat and SSQ)
% marginal estimate (hat and SSQ)
Xd = PX{end}*x;
ssD = trace(Xd*Xd');
res.ExpVar_D_tot = ssD;
res.Xd_crude = Xd;

% projection onto marginal P2
Xd = P2*x;
ssD = trace(Xd*Xd');
res.ExpVar_D_marginal = ssD;
res.Xd_marginal = Xd;


E3f = (eye(n) - P1 - P2)*x;
EtE = E3f'*E3f;  
E2f = (eye(n) - P1)*x;
sse3f = trace(EtE);
diagEtEm = diag(EtE); 

sse2f = trace(E2f'*E2f);

Fm      = ((sse2f - sse3f)      /(df3f - df2f))/(sse3f/(n-df3f));

res.marginalSSQ = max(sse2f - sse3f,0);
res.modelSSE = sse3f;
res.nperm = nperm;
res.rank = df3f; % model rank
res.totVar = trace(x*x');
res.Fmodel = Fm;
res.dftop = (df3f - df2f);
res.dfbottom = n-df3f;

if df3f==df2f
    res.p = 1;
    res.p_F = 1;
    res.p_F_est = NaN;
    return
end


%% setup exchangability kernel probabilty
% predict the class member ship given the propmatch matrix

if prop_matching == 1
    [Pz_x, pYhat, propSWAP] = get_exchangability_kernelprobabilty(propMcov,D{end});
    propSWAP = TSnorm( propSWAP);
end
% permutation
diagEtEperm = []; 
for i=1:nperm
    if isnested==1 && iscrossed==0
        % testing of a systematic effect nested within a random one
        id = nestedrandperm(nestedfac,classD);
    elseif isnested==0 && iscrossed==1
        % make missing within nestedfac
        id = crossedrandperm(nestedfac,classD);
    else
        % Normal ANOVA design _total_ permutation
        %%% propensity matching
        if prop_matching==1
            %%
            
            %rng(123)
            %id0 = permutesampler(Pz_x);
            %rng(123)
            %id = pairwiseswapper(propSWAP);
            %%
            %for jj = 1:100
            id = permutesampler(propSWAP);
            
            %prp(jj,:) = [propSWAP(2,id(2)) mean(propSWAP(2,:))];
            %end
            
            %
            %plot(Pz_x(:),propSWAP(:),'.'); shg
            % compare
            %[id id0]
            %%
            %[u s v ] = svds(auto(auto(Pz_x')'),2);
            %id2 = abs(D{end}(id,1) - D{end}(:,1));
            %scatter(u(:,1),u(:,2),[],id2); shg
            %
            
        else
            %%% all equally likely
            id = randperm(n);
        end
    end
    
%    if unique(D{end}(id,1) - D{end}(:,1))==0
%        'no samples are swapped btw groups - I.e. defacto similar test'
%    end
    
    [ss3fp(i) diagEtEperm(i,:)] = getSSE(P2,P1,defM1,x,id,contvar);
    %ss3fp(i) = getSSE(P2,P1,defM1,x,id2,contvar);
    %ss3fp_pm(i) = getSSE(P2,P1,defM1,x,id1,contvar);
    
    %if contvar==0
    %    P2p = (defM1*P2(id,id))*pinv(defM1*P2(id,id));
    %elseif contvar==1
    %    P2p = P2(id,id)*defM1;
    %end
    %E3fp = (eye(n) - P1 - P2p)*x;
    %ss3fp(i) = trace(E3fp'*E3fp);
    
    Fp(i) = ((sse2f - ss3fp(i))/(df3f - df2f))/(ss3fp(i)/(n-df3f));
    
    % KILL the loop if the p-value is screamingly non-significant
    if sum(Fp>Fm) > 2*alpha_cut*nperm
        length(Fp)
        break
    end
end

if nperm==0
    res.p = NaN;
    res.permutationSSE = NaN;
    res.Fperm = NaN;
    res.p_F = NaN;
    res.diagEtEm = NaN;
    res.diagEtEperm = NaN; 
    res.individual_variable_pv_freq  = NaN; 
else
    indivP = diagEtEperm - ones(nperm,1)*diagEtEm';
    indivP  = (1 + sum(indivP<0))/(1+nperm); 
    res.nperm = length(Fp);
    res.permutationSSE = ss3fp;
    res.diagEtEm = diagEtEm;
    res.diagEtEperm = diagEtEperm;
    res.individual_variable_pv_freq = indivP;   
    res.p = (sum(ss3fp<=sse3f)+1)/(res.nperm+1);
    res.Fperm = Fp;
    res.p_F = (sum(Fp>Fm)+1)/(res.nperm+1);
end

% estimate F distr parameters
try
    %
    vv = mle(Fp,'pdf',@fpdf,'start',[res.dftop res.dfbottom]);
    res.dftop_est = vv(1);
    res.dfbottom_est = vv(2);
    res.p_F_est = 1 - fcdf(Fm,vv(1),vv(2));
    %vv = mle(Fp,'pdf',@ncfpdf,'start',[1 3 0]);
catch
    res.dftop_est = NaN;
    res.dfbottom_est = NaN;
    res.p_F_est = NaN;
end

res.prop_matching = prop_matching;
%res.sse_perm =  ss3fp;
%res.sse_perm_propmatch =  ss3fp_pm;
% res.timeconsumption = tm;

%
[u s v0] = svds(res.Xd_crude,1);

% inverse 
z = abs(norminv(res.individual_variable_pv_freq/2)); 
[u s v1] = svds(res.Xd_crude*diag(z),1);
res.v0 = v0; 
res.v1 = v1;  
%plot(v1,v0,'*'); shg
%
%plot(res.individual_variable_pv_freq,z,'.'); shg


if plotit==1
    plotANOVApermutationtest(res)
end

function [E, sse] = calcsse(D,x)

if nargin==1 % only x;
    E = D;
else
    E = x - D*pinv(D)*x;
end

[n p] = size(E);

if n<p; sse= trace(E*E');
else
    sse = trace(E'*E);
end



function [E, sse] = calcsse2(DefM,x)

if nargin==1 % only x;
    E = DefM;
else
    E = DefM*x;
end

[n p] = size(E);

if n<p; sse= trace(E*E');
else
    sse = trace(E'*E);
end



function Xn = normaliz(X)
n = size(X,2);
Xn  = X./(sum(X,2)*ones(1,n));


function [sse diagEtE] = getSSE(P2,P1,defM1,x,id,contvar)

n = size(P2,1);
if contvar==0
    P2p = (defM1*P2(id,id))*pinv(defM1*P2(id,id));
elseif contvar==1
    P2p = P2(id,id)*defM1;
end
E3fp = (eye(n) - P1 - P2p)*x;
EtE = E3fp'*E3fp;
diagEtE = diag(EtE ); 
sse = trace(EtE);
