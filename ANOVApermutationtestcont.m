function res = ANOVApermutationtestcont(x,PX,D,nperm,contvar,prop_matching,slopecorr,classD,nestedfac)

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
if nargin<8
    % Normal ANOVA type model, NO mixed effects
    isnested = 0;
    iscrossed = 0;
else
    d1 = mkdesignmatrix(classD);
    d2 = mkdesignmatrix(nestedfac);
    isnested = (rank(d2)==rank([d1 d2]))+0;
    
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

n = size(x,1);

%Xhat_crude = PX{end}*x;
[P1 P2] = blockwisehat2(D(1:end-1));
%DMorth = (eye(n) - P1 - P2)*PX{end};     % orthogonalize Kernel with _all others_
%Xhat_marginal= DMorth*(eye(n) - P1 - P2)*x;  % Hat on residuals after removal of _all others_
%a = pinv(Xhat_marginal(:))*Xhat_crude(:);        % make slope correction in reference to crude estimate
%Xhat_marginal = Xhat_marginal*a;
[Xhat_crude,Xhat_marginal,DMorth,a] = projectContKernel(x,PX{end},P1+P2,slopecorr);

df3f = trace(DMorth*a + P1 + P2);
df2f = trace(P1 + P2);

% crude estimate (hat and SSQ)
% marginal estimate (hat and SSQ)
Xd = Xhat_crude;
ssD = trace(Xd*Xd');
res.ExpVar_D_tot = ssD;
res.Xd_crude = Xd;

% projection onto marginal P2
Xd = Xhat_marginal;
ssD = trace(Xd*Xd');
res.ExpVar_D_marginal = ssD;
res.Xd_marginal = Xd;

E3f = (eye(n) - P1 - P2)*x - Xhat_marginal;
E2f = (eye(n) - P1 - P2)*x;

sse3f = trace(E3f'*E3f);
sse2f = trace(E2f'*E2f);

Fm      = ((sse2f - sse3f) /(df3f - df2f))/(sse3f/(n-df3f));

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


% setup exchangability kernel probabilty
% predict the class member ship of the last design in D given all others
%Y = D{end};
%Yhat = zeros(size(Y));
%DD = ones(n,1);
%for i=1:(length(D)-1)
%    DD = [DD D{i}];
%end

%for i=1:size(Y,2)
%    b = glmfit(DD,Y(:,i),'binomial','link','logit');
%    Yhat(:,i) = glmval(b,DD,'logit') + 0.01;
%end

%Yhat = normaliz(Yhat);

% between all pairs, calculate the kullback-leiber divergence metric.
%for j=1:n
%    KLdist(:,j)=KLDiv(Yhat,Yhat(j,:));
%end
% transform the metric into probabilities of exchangegability

% KLD is asymptotically distributed as a scaled (non-central) chi-square
% with one degree of freedom (ref: British Journal of Mathematical and Statistical Psychology (2011), 64, 291?309)
%hist(KLdist(:),100);
%
%[~,idd] =  sort(Yhat(:,1));
%Pz_x = 1 - chi2cdf(KLdist,1);
%Pz_x = normaliz(Pz_x);
%Pz_x = Pz_x(idd,idd);

%for ii=1:10000
%    ID(ii,:) = permutesampler(Pz_x);
%end


%% sum(Pz_x ,2)
%plot(KLdist(:),Pz_x(:),'*'); shg
%

%gm = exp(mean(log(Yhat),2));
%Yhat = diag(1./gm)*Yhat;
%Dst = squareform(pdist(Yhat));

%Pz_x = 1 - normcdf(Dst / (sqrt(2/n)));
%Pz_x = normaliz(Pz_x)
%

% permutation
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
            id = permutesampler(Pz_x);
        else
            %%% all equally likely
            id = randperm(n);
        end
    end
    
    
    
    [~,Xhat_marginalP] = projectContKernel(x,PX{end}(id,id),P1+P2, slopecorr);
    E3f = (eye(n) - P1 - P2)*x - Xhat_marginalP;
    ss3fp(i) = trace(E3f'*E3f);
    
    Fp(i) = ((sse2f - ss3fp(i))/(df3f - df2f))/(ss3fp(i)/(n-df3f));
    
    % KILL the loop if the p-value is screamingly non-significant
    if sum(Fp>Fm) > 2*alpha_cut*nperm
        length(Fp)
        break
    end
end

res.S2Ndataspace = NaN; % for now just leave this out. 
res.S2Nsca = NaN; 

if nperm==0
    res.p = NaN;
    res.permutationSSE = NaN;
    res.Fperm = NaN;
    res.p_F = NaN;
    res.p_sca = NaN; 
else
    res.nperm = length(Fp);
    res.permutationSSE = ss3fp;
    res.p = (sum(ss3fp<sse3f)+1)/(res.nperm+1);
    res.Fperm = Fp;
    res.p_F = (sum(Fp>Fm)+1)/(res.nperm+1);
    res.p_sca = NaN; % for now just leave this out. 
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


function sse = getSSE(P2,P1,defM1,x,id,contvar)

n = size(P2,1);
if contvar==0
    P2p = (defM1*P2(id,id))*pinv(defM1*P2(id,id));
elseif contvar==1
    P2p = P2(id,id)*defM1;
end
E3fp = (eye(n) - P1 - P2p)*x;
sse = trace(E3fp'*E3fp);
