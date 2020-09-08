function [P1 P2] = blockwisehat2(D,contvar)
%
% THE BLOCK update of full model hat matrix from a sequence of crude
% hat matrices.
%%% input D: list of DESIGN MATRIXES
%%% output: two ortogonal projection matrices: P2 the marginal of the last
%%% matrix, P1 = all the other effects in the model

nterms = length(D);
n = size(D{1},1);
P1 = ones(n,1)*pinv(ones(n,1)); 

e = 1e-9;
for i=1:nterms
      DefM1 = eye(n) - P1;
      BorthP1 = DefM1*D{i};
      if norm(BorthP1)>e
      P2 = (BorthP1)*pinv(BorthP1);
      else
          P2 = zeros(n);
      end
      
    % check if it is a projection matrix or not
    %isprjm = checkProjM(P2,1e-6);
    
  
    %PX2 = W2;
    % update P2 (marginal effect given P1)
    %if isprjm
    %    P2 = (DefM1*P2)*pinv(DefM1*P2);
    %else
        %P2 = DefM1*P2*DefM1;
    %    P2 = P2*DefM1; 
%        P2 = DefM1*P2; 
    %end
    %PX1 = (DefM2*W1)*pinv(DefM2*W1);
    % update P1
    P1 = P1 + P2;
end
P1 = P1 - P2;


function flag = checkProjM(P2,tol)

sumabsP2 = sum(abs(P2(:))); 
D = P2*P2 - P2;
sumabs = sum(abs(D(:)));
flag = (sumabs / sumabsP2)  <tol;
