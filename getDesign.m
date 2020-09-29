function D = getDesign(k,rep,thr)

% make nested mixed design with null fixed factor but significant random factor
%k = 30; % individuals
%rep = 5; % repeated measures
D1 = fullfact([rep k]);
D2 = (D1(:,2)<=(k/2))+1;
D = [D1(:,2) D2 D1(:,1)]; % ID, condition, trt
%lbD = {'ID','condition', 'trt'};

% make unbalanced
icout1 = D2==1 & ismember(D1(:,1),[1 2]) & rand(length(D2),1)>thr;
icout2 = D2==2 & ismember(D1(:,1),[4 5]) & rand(length(D2),1)>thr;
icout = icout1 | icout2;
D = D(~icout,:); 
 