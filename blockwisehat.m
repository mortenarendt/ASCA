function [P1 P2] = blockwisehat(PX)
%
% THE BLOCK update of full model hat matrix from a sequence of crude
% hat matrices.
%%% input PX: list of Projection matrices
%%% output: two ortogonal projection matrices: P2 the marginal of the last
%%% matrix, P1 = all the other effects in the model

nterms = length(PX);
%P1 = PX{1};
n = size(PX{1},1);
%if nterms==1
%P2 = zeros(n); 
P1 = ones(n,1)*pinv(ones(n,1)); 
%    P2 = PX{1};
%    return
%end

for i=1:nterms
    P2 = PX{i};
    % check if it is a projection matrix or not
    isprjm = checkProjM(P2,1e-6);
    
    DefM1 = eye(n) - P1;
    %PX2 = W2;
    % update P2 (marginal effect given P1)
    if isprjm
        P2 = (DefM1*P2)*pinv(DefM1*P2);
    else
        %P2 = DefM1*P2*DefM1;
        P2 = P2*DefM1; 
%        P2 = DefM1*P2; 
    end
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
