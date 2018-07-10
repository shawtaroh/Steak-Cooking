%
% chi.m - returns flory-huggins interaction parameter
%
%
function y = chi(P, T, phi)
    y = chi_p( P,T) - (chi_p(P, T)-P.chi_0*ones(size(phi))).*((ones(size(phi))-phi).^2);
    %y=.8
end


    function x = chi_p(P, T)
        exp_term = P.A*exp(P.gamma*(P.T_0*(ones(size(T))+T*P.nu/P.alpha)-P.T_e*ones(size(T))));
       x =  P.chi_pn*ones(size(T)) + (P.chi_pd-P.chi_pn)*ones(size(T))./(ones(size(T))+exp_term);
    end

