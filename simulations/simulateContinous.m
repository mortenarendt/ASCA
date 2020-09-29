close all hidden; clear; clc

clear
p = 20;
% make nested mixed design with null fixed factor but significant random factor
k = 30; % individuals
rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<16)+1;
covariate = randn(size(D1,1),1);
D = [D1(:,2) D2 D1(:,1) covariate]; % ID, condition, trt
lbD = {'ID','condition', 'trt','covar'};

% make unbalanced
% remove 50% of samples D2==1 and D1 in 1 and 2
% remove 50% of samples D2==2 and D1 in 4 and 5
icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>0.5;
icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>0.5;
icout = icout1 | icout2;
sum(icout)
D = D(~icout,:);

M = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1;  0 1 1 0 ];
s = [0.1 0.1 0.1 0.03 0 .2];

%s = [0 0 0 0 0 1];


plotit=0;

opt = ASCAcat('options');
opt.contvar = 4;
opt.contvaraslinear = 1;
opt.nperm = 1000;
opt.showtable = 0;
clear pv;

EXPVAR = nan(50,4); 
%
for j=1%:10
j
% new covariate
D(:,4) = randn(size(D,1),1);
% simulate data 
[X DM Xm] = simulateData(D,M,s,p,opt.contvar);

opt.contvaraslinear = 1;
opt.slopecorr = 0;
reslin = ASCAcat(X,D,M, opt);

opt.contvaraslinear = 0;
opt.slopecorr = 0;
ressm = ASCAcat(X,D,M, opt);

opt.slopecorr = 1;
ressmslc = ASCAcat(X,D,M, opt);

EXPVAR(j,:) = [trace(Xm{4}'*Xm{4}) reslin.Effects{4}.ExpVar_D_marginal ressm.Effects{4}.ExpVar_D_marginal ressmslc.Effects{4}.ExpVar_D_marginal];
end


% % % 
% % % n = size(X,1);
% % % Dlin = D(:,4)*pinv(D(:,4));
% % % 
% % % 
% % % Xhat{1}= DM{4}*X;
% % % Xhattot{1} = Xhat{1};
% % % Xhat{2}= Dlin*X;
% % % Xhattot{2} = Xhat{2};
% % % %%% get rid of the inteferents
% % % [P1 P2] = blockwisehat2(DM([1 3 5]));
% % % Xhat{3}= DM{4}*(eye(n) - P1 - P2)*X;
% % % Xhattot{3} = (P1 + P2)*X + Xhat{3};
% % % % correlate crude with marginal to update marginal with a slope
% % % a = pinv(Xhat{3}(:))*Xhat{1}(:);
% % % Xhat{7} = Xhat{3}*a;
% % % Xhattot{7} = (P1 + P2)*X + Xhat{7};
% % % 
% % % idd = randperm(n);
% % % DM4orth = (eye(n) - P1 - P2)*DM{4};     % orthogonalize Kernel with _all others_
% % % Dlinorth = (eye(n) - P1 - P2)*D(:,4);   %
% % % Dlinorth = Dlinorth *pinv(Dlinorth );
% % % Xhat{5}= DM4orth*(eye(n) - P1 - P2)*X;  % Hat on residuals after removal of _all others_
% % % Xhattot{5} = (P1 + P2)*X + Xhat{5};     % The full hat
% % % 
% % % %Xh = (eye(n) - P1 - P2)*DM{4}(idd,:)*(eye(n) - P1 - P2)*X;
% % % %dm = (eye(n) - P1 - P2)*DM{4}*(eye(n) - P1 - P2);
% % % %Xh2 = dm(idd,idd)*X;
% % % %plot(Xh2(:), Xhat{5}(:),'.'); shg
% % % %%
% % % a = pinv(Xhat{5}(:))*Xhat{1}(:);        % make slope correction in reference to crude estimate
% % % Xhat{8} = Xhat{5}*a;
% % % a
% % % Xhattot{8} = (P1 + P2)*X + Xhat{8};
% % % 
% % % 
% % % Xhat{6}= Dlinorth*(eye(n) - P1 - P2)*X;
% % % Xhattot{6} = (P1 + P2)*X + Xhat{6};
% % % Xhat{4} = Dlin*(eye(n) - P1 - P2)*X;
% % % Xhattot{4} = (P1 + P2)*X + Xhat{4};
% % % 
% % % 
% % % 
% % % 
% % % 
% % % for i=1:8
% % %     %E = Xhat{i} - Xm{4};
% % %     E = X - Xhattot{i};
% % %     sse(i) = trace(E'*E);
% % %     if plotit==1
% % %         subplot(3,3,i);
% % %         %plot(Xm{4}(:),Xhat{i}(:),'.');
% % %         plot(X(:) - Xm{end}(:),Xhattot{i}(:),'.');
% % %         hold on; abline(1,0); shg
% % %         title(num2str(sse(i),4))
% % %         title([LB{i} ' see = ' num2str(sse(i),4)])
% % %     end
% % % end
% % % SSE(j,:) = sse;
% % % %end
% % % 
% % % 
