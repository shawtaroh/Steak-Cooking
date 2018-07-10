%
% xi - returns local aspect ratio
%
function y = xi(P,T)
%    y = 1 + .1/(1+P.A*exp(P.gamma*((T*(P.T_D-P.T_0)+P.T_0)-P.T_e*ones(size(T,1),size(T,2)))))
        exp_term = P.A*exp(P.gamma*(P.T_0*(ones(size(T))+T*P.nu/P.alpha)-P.T_e*ones(size(T))));
       y =  ones(size(T)) + (.25)*ones(size(T))./(ones(size(T))+exp_term);
end