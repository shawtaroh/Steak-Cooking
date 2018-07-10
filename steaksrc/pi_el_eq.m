%
%   pi_el.m - returns elastic pressure for variable phi_0
%
%
function y = pi_el_eq(P,T,phi)
        y = P.beta_el*(ones(size(T))+P.nu/P.alpha*T).*(phi/2);
end
