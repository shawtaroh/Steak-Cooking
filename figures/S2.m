
    RHS_1 = P.omega^2*ones(size(S2.phiold))./(P.omega*(ones(size(S2.phiold))-S2.phiold)+ones(size(S2.phiold))).^2;
    RHS_1 = RHS_1.*S2.phi_cdx.*S2.T_cdx;
    RHS_2 = P.omega*ones(size(S2.phiold))./(P.omega*(ones(size(S.phiold))-S2.phiold)+ones(size(S2.phiold))).*S2.T_xx;
    
    RHS_3 = (P.omega-1)*(S2.phi_cdy.*S2.T_cdy);
    RHS_4 = (ones(size(S2.phiold))-(1-P.omega)*S2.phiold).*S2.T_yy;
    
    RHS_5 = (S2.phiold-ones(size(S2.phiold))).*(S2.w(:,:,1).*S2.T_cdx+S2.w(:,:,2).*S2.T_cdy);
    
    RHS_6 = -(P.alpha*ones(size(S2.T))+P.nu*S2.T).*(S2.phi-S2.phiold)/P.dt;
    