function results = ASCAcat(X,dep,M,deplb,options)
%
% results = ASCAcat(X,dep,M,deplb,options)
% X (n,p) prepreocessed input data
% dep (n,m) matrix or vector with design factors
% M (k,m) a 0/1 model matrix indicating which of the m-factors that are included in each of the k model terms.
% E.g. M = [1 0; 0 1; 1 1] refers to a design with two factors and a model
% deplb (m,1) labels for the individual design vectors
% with the two main effects (row 1 and 2) and interaction between those
% (row 3)
% options = ASCAcat('options') provides a list of algorithmic options.

%%%%% Check the input %%%%%%

if nargin==1 && ischar(X)
    % make default option structure
    results.nperm = 1000; % number of permutations
    results.deflate = 1; % deflate to get marginal
    results.prop_matching = 0; % do not do prop_matching in the randomization
    results.rndfac = []; % Default is no mixed effects.
    results.contvar = []; % Defalt is all design variables being factors
    results.showtable = 1; % Print ANOVA table to screen
    results.slopecorr = 1; % make a scalar slope correction in the continous kernel estimation
    results.contvaraslinear = 0;
    return
end

nf = size(dep,2);
% check if options are provided
if nargin==2
    %%% results = ASCAcat(X,dep)
    options = ASCAcat('options');
    M = eye(nf);
    deplb = getfactorlabels(nf);
end


if nargin==3
    
    if isstruct(M)
        %%% results = ASCAcat(X,dep,options)
        options = M;
        M = eye(nf);
        deplb = getfactorlabels(nf);
    else
        %%% results = ASCAcat(X,dep,M)
        options = ASCAcat('options');
        deplb = getfactorlabels(nf);
    end
end


if nargin==4
    %%% results = ASCAcat(X,dep,M, options)
    if isstruct(deplb)
        options = deplb;
        deplb = getfactorlabels(nf);
    else
        %%% results = ASCAcat(X,dep,M, deplb, options)
        options = ASCAcat('options');
    end
end

% expand nperm to vector of length == number of terms
if numel(options.nperm)==1
    nterms = size(M,1);
    options.nperm = ones(1,nterms)*options.nperm;
end



% remove samples with missing design info
if size(dep,2)==1
    ic = ~isnan(dep);
else
    ic = sum(isnan(dep'))==0;
end
nout = sum(ic==0);
if nout>0
    disp(['#' num2str(nout) ' Samples have missing info in design and are removed'])
end
X = X(ic,:);
dep = dep(ic,:);

n = size(X,1);
nfac = size(dep,2);

if nfac~=size(M,2)
    error(['The design (M) needs to include all factors (#factors = ' num2str(nfac)]);
end

if isdataset(X)
    X = X.data;
end


%%%% START the modelling here %%%%%%%

% based on which factors that are considered random, find which terms of
% the model that are random:
rndterms = find(sum(M(:,options.rndfac),2)>0);

results.model = M;
results.random_factors = options.rndfac;
results.random_terms = rndterms;
results.continous_factors = options.contvar;
results.data = X;
results.dependent_var = dep;
results.dependent_var_lb = deplb;
results.prop_matching = options.prop_matching;
results.options = options;

% setup design Matrices D (dummy for categorical or kernel weighting for covariates)
% and deflation matrices for each term in the model
nterms = size(M,1);
order = sum(M');
contterm = zeros(1,nterms);
for i=1:nterms
    nf = sum(M(i,:));
    iddep = find(M(i,:));
    if nf==1 % Main effects;
        if ismember(iddep,options.contvar)
            % continous variable
            if options.contvaraslinear ==0
                PX{i} = setCovariateKernel(dep(:,iddep));
                D{i} = dep(:,iddep);
            else
                PX{i} = [ones(n,1) dep(:,iddep)]*pinv([ones(n,1) dep(:,iddep)]);
                D{i} = [ones(n,1) dep(:,iddep)];
            end
            contterm(i) = 1;
            
        else
            % categorical variable
            D{i} = mkdesignmatrix(dep(:,iddep));
            contterm(i) = 0;
            PX{i} = D{i}*pinv(D{i});
        end
        
        IDvec(:,i) = dep(:,iddep);
        lbD(i) = deplb(iddep);
    elseif nf==2 % two-way interaction effects.
        [D{i},~,~,IDvec(:,i)] = mkdesignmatrix(dep(:,iddep(1)),dep(:,iddep(2)));
        lbD(i) = cellstr([char(deplb(iddep(1))) ' x ' char(deplb(iddep(2)))]);
        PX{i} = D{i}*pinv(D{i});
    elseif nf>2
        [D{i},IDvec(:,i)] = mkdesignmatrixHigherOrder(dep(:,iddep));
        PX{i} = D{i}*pinv(D{i});
        l = [];
        for ii = 1:length(iddep)
            if ii==length(iddep)
                l = [l char(deplb(iddep(ii)))];
            else
                l = [l char(deplb(iddep(ii))) ' x '];
            end
        end
        lbD(i) = cellstr(l);
    end
end
%
% find nested structure and estimate marginal design matrices
m = FindNestedStructure(M);
m2 = FindNestedStructure_DoE(PX,M);
% update the nested structure for the testing
m = (m + m2)>0;

for i=1:nterms
    %if i==6
    %    whos
    %end
    % nested structure
    %%% mi := which terms are coarser than term i
    mi = m(i,:);
    mi(i) = 0;
    %%% mj := which terms include term i
    mj = m(:,i);
    mj(i) = 0;
    nestedeff = find(mi==1);
    
    DD = [];
    % if directly nested under factors (e.g interactions under maineffect
    % or coarser main effect)
    for ii=1:length(nestedeff)
        DD = [DD D{nestedeff(ii)}];
    end
    PM = [];
    DefM = eye(n);
    if options.deflate == 1
        % if any covariates that are correlated (for unbalanced designs)
        whichother = find(mj==0);
        whichother = whichother(~ismember(whichother ,i));
        
        for ii=1:length(whichother)
            DD = [DD D{whichother(ii)}];
            PM = [PM D{whichother(ii)}];
            % covariate
            if ismember(whichother(ii) ,options.contvar) & options.contvaraslinear==0
                DefM = DefM*(eye(n) - D{whichother(ii)});
            else
                DefM = DefM*(eye(n) - D{whichother(ii)}*pinv(D{whichother(ii)}));
            end
            
        end
        
    end
    
    % Calculate Marginal effects.
    if contterm(i)==1
        Dmarg{i} = D{i};
    else
        Dmarg{i} = DefM*D{i};
        Dmarg{i} = (eye(n) - PM*pinv(PM))*D{i};
    end
    df(i) = rank(Dmarg{i}); % marginal degrees of freedom
    
    whichrndterm = intersect(rndterms,find(mj));
    
    if ismember(rndterms,i)
        % if term is random, then DO NOT test it, but still return its SS values
        res{i} = ANOVApermutationtest(X,...
            PX(i),D(i),0,0,0);
        res{i}.israndomfactor = 1;
    else
        % if a term is nested within a random factor (nested design), use this
        % factor in the permutation testing
        if ~isempty(whichrndterm)
            % Which coarser factor? whichrndterm
            % get the indicator vector for this term
            %res = ANOVApermutationtest(x,PX,D,nperm,contvar,prop_matching,classD,nestedfac)
            res{i} = ANOVApermutationtest(X,...
                PX([unique([nestedeff whichother(:)']) i]),...
                D([unique([nestedeff whichother(:)']) i]),...
                options.nperm(i),0,0,IDvec(:,i),IDvec(:,whichrndterm)); % calculate statistical inference.
        elseif ~isempty(options.rndfac) & isempty(whichrndterm)
            % if there is mixed effects but crossed with systematic factors, use this in the permutation testing
            res{i} = ANOVApermutationtest(X,...
                PX([unique([nestedeff whichother(:)']) i]),...
                D([unique([nestedeff whichother(:)']) i]),...
                options.nperm(i),0,0,IDvec(:,i),IDvec(:,rndterms)); % calculate statistical inference.
        else
            if contterm(i)==1
                if options.contvaraslinear==0
                    res{i} = ANOVApermutationtestcont(X,...
                        PX([unique([nestedeff whichother(:)']) i]),...
                        D([unique([nestedeff whichother(:)']) i]),...
                        options.nperm(i),contterm(i), options.prop_matching, options.slopecorr); % calculate statistical inference.
                else
                    res{i} = ANOVApermutationtest(X,...
                        PX([unique([nestedeff whichother(:)']) i]),...
                        D([unique([nestedeff whichother(:)']) i]),...
                        options.nperm(i),contterm(i), options.prop_matching); % calculate statistical inference.
                end
                res{i}.israndomfactor = 0;
            else
                % normal vanilla ANOVA stats
                res{i} = ANOVApermutationtest(X,...
                    PX([unique([nestedeff whichother(:)']) i]),...
                    D([unique([nestedeff whichother(:)']) i]),...
                    options.nperm(i),contterm(i), options.prop_matching); % calculate statistical inference.
                res{i}.israndomfactor = 0;
            end
            
        end
    end
    %Xmarginal{i} = res{i}.Xd_marginal;
    % Store the results.
    SS(i)       = res{i}.marginalSSQ;
    SScrude(i)  = res{i}.ExpVar_D_tot;
    p(i) = res{i}.p_F;
    pest(i) = res{i}.p_F_est;
end

for i=1:nterms
    idd = true(nterms,1);
    idd(i) = false;
    idd = find(idd);
    DD = [];
    for ii=1:length(idd)
        DD = [DD Dmarg{idd(ii)}];
    end
    if ~isempty(DD)
        Dmarg2{i} = (eye(n) - DD*pinv(DD))*D{i};
    end
end

%%
DD = []; id = [];
for i=1:nterms
    %     DD = [DD Dmarg{i}];
    %     DD =[DD Dmarg2{i}];
    [u s v] = svds(Dmarg{i}, rank(Dmarg{i}));
    idkeep = diag(s)>1e-7;
    DD = [DD u(:,idkeep)];
    %     id = [id ones(1,size(D{i},2))*i];
    id = [id ones(1,sum(idkeep))*i];
end

%%
SStot = trace(X'*X);
[P1 P2] = blockwisehat(PX);
Xhat = (P1 + P2)*X;
E = X - Xhat;
SSE = trace(E'*E);

% check balancedness quality of each marginal effect
for i=1:length(Dmarg)
    practicalrank = rank(Dmarg{i});
    % can we estimate it?
    ss = svds(Dmarg{i}'*Dmarg{i},practicalrank);
    %if i==3
        %[i rank(D{i}) practicalrank ss']
    %end
    balancedness(i) =  ss(end) / ss(1);
    res{i}.balancedness = balancedness(i);
    res{i}.balancedness_singval = ss';
end

results.balancedness = balancedness;

% % Check SS consistency
% sum(ssdd) + SSE, SStot
% correct ssdd such that sum(ssdd) = SStot - SSE

%%ssdd = ssdd * ((SStot - SSE) / sum(ssdd));
ssdd = SS;
SS = [ssdd SSE SStot ]';
SSmarginal = [ssdd/SStot SSE/SStot];
SSmarginal = [SSmarginal sum(SSmarginal )]';

sscrude = [SScrude/SStot SSE/SStot];
sscrude = [sscrude sum(sscrude)]';

lbSS = [lbD {'Error'} {'Total'}]';
%cellstr(num2str(SS / SS(end),'%1.3f')) ...

atab = [lbSS cellstr(num2str(SS,3)) ...
    cellstr(num2str(sscrude,'%1.3f')) ...
    cellstr(num2str(SSmarginal,'%1.3f')) ...
    [cellstr(num2str(df',3)); {' '};{' '}] ...
    [cellstr(num2str(balancedness',3)); {' '};{' '}] ...
    [cellstr(num2str(p',3)); {' '};{' '}] ...
    [cellstr(num2str(pest','%1.3g')); {' '};{' '}]];
atab = [{'Factor' 'SS' 'SS/SStot' 'SSmarg/SStot' 'DF' 'BalNess' 'P-value-freq' 'P-value-Fest'};atab];

results.Effects = res;
results.ANOVAtab = atab;
%results.Xmarginal = Xmarginal;
results.E = E;
%results.D = DD;
%results.reg = b;
%results.regcoef_termID = id;


if options.showtable==1
    try
        figh=disptable(results.ANOVAtab,'ASCA',['ANOVA table (#perm ' num2str(options.nperm) ')']);
        results.ANOVAtab_handle = figh;
    end
end

% plot results
% plotASCAresultsV2(results);
results = makeScores(results);

function m = FindNestedStructure(M)
% takes model matrix M (terms by design variables) and returns a terms by
% terms indicating which terms that are nested for any specific factor.
% The main test:  no overlap between terms included in t(i) and not in t(j)
m1 = (M*M'); %> 0
nfac = (M*ones(size(M,2),size(M,1)))';
m = m1>=nfac;

function m = FindNestedStructure_DoE(D,M)
% given the specific design (reflected by D) and the model in M, figure out
% which terms that are a finner version of any other.

nterms = size(M,1);
m = ones(nterms);
for i=1:nterms
    %for the i'the term, check if it is nested within any of the other of terms
    idnoti = true(1,nterms);
    idnoti(i) = false;
    for j = find(idnoti)
        %Dj = D{j};
        %m(i,j) = rank([D{i} D{j}])==rank(D{i})
        m(i,j) = checkmapping(D{j},D{i},1e-6);
    end
    
end

%%

function flag = checkmapping(D1,D5, tol)

%D5 = D{5}; D1 = D{1};
n = size(D5,1);
% can D1 be created as a linear mapping of D5?
% D1 = D5*B
E = (eye(n) - D5*pinv(D5))*D1;
flag = trace(E'*E)<tol;

function [D, cl] = mkdesignmatrixHigherOrder(dep)

nterms = size(dep,2);

% make first two order interaction
[D, ~,~, cl1] = mkdesignmatrix(dep(:,1),dep(:,2));

if nterms>2
    for i = 3:nterms
        [D, ~,~,cl] = mkdesignmatrix(cl1,dep(:,i));
    end
end


function deplb = getfactorlabels(nf)
deplb = cellstr([repmat('Factor',nf,1) num2str((1:nf)')])';

