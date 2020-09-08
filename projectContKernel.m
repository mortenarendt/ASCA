function [Xhat_crude,Xhat_marginal,DMorth,a] = projectContKernel(x,PX,PP,slopecorr)

n = size(x,1);
Xhat_crude = PX*x;
DMorth = (eye(n) - PP)*PX;     % orthogonalize Kernel with _all others_
Xhat_marginal= DMorth*(eye(n) - PP)*x;  % Hat on residuals after removal of _all others_
if slopecorr==1
    a = pinv(Xhat_marginal(:))*Xhat_crude(:);        % make slope correction in reference to crude estimate
else
    a = 1;
end
Xhat_marginal = Xhat_marginal*a;


%Xhat_crude = PX{end}*x;
%[P1 P2] = blockwisehat2(D(end-1));
%DMorth = (eye(n) - P1 - P2)*PX{end};     % orthogonalize Kernel with _all others_
%Xhat_marginal= DMorth*(eye(n) - P1 - P2)*x;  % Hat on residuals after removal of _all others_
%a = pinv(Xhat_marginal(:))*Xhat_crude(:);        % make slope correction in reference to crude estimate



