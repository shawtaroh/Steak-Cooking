% 
% j_evap.m - returns non-dim term 8.2.2 
% 
function y = j_evap(P,T,phi,t)

    c_0 = P.c_0;
    if(P.T_D == 225+273)
        c_0 = c_0_225(P,t);
    end
    if(P.T_D == 175+273)
        c_0 = c_0_175(P,t);
    end
    
    y = P.lambda_bar/P.c_0*c_0*((P.p_sat_0*P.M_f*P.X_m)/P.R/c_0*ones(size(T))./((P.T_D-P.T_0)*T+P.T_0*ones(size(T))).*(1/P.X_m*ones(size(T))-P.rho_s/P.rho_f*phi./(ones(size(T))-phi))); 
    y = y.*exp(17.27*((P.T_D-P.T_0)*(T)+(P.T_0-273)*ones(size(T)))./((P.T_D-P.T_0)*T+(P.T_0-35.86)*ones(size(T)))); 
    y = y-P.lambda_bar/P.c_0*c_0*ones(size(T)); 
    y = zeros(size(T));
end