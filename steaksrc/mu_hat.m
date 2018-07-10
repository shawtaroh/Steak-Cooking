%
%   mu.m - Returns mu_hat(T)
%
%
function y = mu_hat(P,T)

y = mu(T*(P.T_D-P.T_0)+P.T_0*ones(size(T)));
%y = 2.414e-5 * 10^(247.8*ones(size(T))./(T-140*ones(size(T))));
y = y/P.mu_0;

end