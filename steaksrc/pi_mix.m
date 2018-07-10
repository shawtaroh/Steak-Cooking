%
%   pi_mix.m - returns osmotic pressure
%
%
function y = pi_mix(P,T,phi)
    y = P.beta_mix*(ones(size(T))+P.nu/P.alpha*T).*(log(1-phi)+phi+chi(P,T,phi).*(phi.^2));
end