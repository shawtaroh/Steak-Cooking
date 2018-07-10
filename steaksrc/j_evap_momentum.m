% 
% j_evap_momentum.m - returns non-dim term 
% 
function y = j_evap_momentum(P,T,phi,t) 


    c_0 = P.c_0;
    if(P.T_D == 225+273)
        c_0 = c_0_225(P,t);
    end
    if(P.T_D == 175+273)
        c_0 = c_0_175(P,t);
    end
    
    A = P.t_0/P.l/P.rho_f*c_0*P.h*((P.D_a/P.k_a)^2/P.rho_a/P.c_a)^(1/3);
    y = A*((P.p_sat_0*P.M_f*P.X_m)/c_0*ones(size(T))./((P.T_D-P.T_0)*T+P.T_0*ones(size(T))).*(1/P.X_m*ones(size(T))-P.rho_s/P.rho_f*phi./(ones(size(T))-phi))); 
    y = y.*exp(17.27*((P.T_D-P.T_0)*(T)+(P.T_0-273)*ones(size(T)))./((P.T_D-P.T_0)*T+(P.T_0-35.86)*ones(size(T)))); 
    y = y-A*ones(size(T)); 
    y=zeros(size(T));
end