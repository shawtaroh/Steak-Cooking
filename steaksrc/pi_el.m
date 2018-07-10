%
%   pi_el.m - returns elastic pressure
%
%
function y = pi_el(P,T,phi)
    %y = P.beta_el_vm*(1+P.nu\P.alpha*T)*(phi^(1/3)*P.phi_0_vm^(2/3)-phi/2);
        y = P.beta_el*(ones(size(T))+P.nu/P.alpha*T).*(phi.^(1/3)*P.phi_0^(2/3)-phi/2);

end
