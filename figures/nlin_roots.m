%
% nlin_roots.m - solves pi_sw=0
%
P = setdefaultparams_steak();

%S.T*(P.T_D-P.T_0)+ones(size(S.T))*(P.T_0-273)
T = 0*ones(size(S.phi(end/2,:)));
phi = S.phi(end/2,:);
pi_sw = P.mu_0/P.t_0*(pi_el(P,T,phi) + pi_mix(P,T,phi))
% phi_0 = 20.5, phi_final = .54
P.chi_pd=.9;
P.chi_pn=.7;

%P.chi_pd=.92;
%P.chi_pn=.72;


%P.phi_0 = 0.174029030084110;
% for phi*(T=0) = %17.403, phi*(T=1)=%43.472
step=1e-2
T = 0:step:1;
x0 = .4:step:.8;
for i =1:length(T)
    count =1;
    for j = 1:length(x0)
        [x,fval,exitflag,output]=root(P,T(i),x0(j));
        if(exitflag>0 && imag(x(1)) ==0 && imag(fval(1)) ==0)
        phi_root(count) = x(1);
        T_intt(count) = T(i);
        x_g(count) = x0(j);
        ffval(count,:) = fval;
        count = count+1;
        end
    end
    phi{i} = phi_root;
    T_gg{i} = x_g;
    fffval{i} = ffval;
    T_int{i} = T_intt;
    clear ('T_intt','T_root','phi_root')
end


for k = 1:length(phi)
    x(k)=median(T_int{k});
    y(k)=median(phi{k});
end


return
plot(T,phi,'LineWidth',2)
title('\phi Equilibrium');

return;
P.chi_pd=1;
P.chi_pn=.05;
P.chi_0=.5;
T = 273:.1:373;
for i =1:length(T)
    c(i) = chi_p(P,T(i));
    phi(i) = root_dim(P,T(i));
end
plot(T,phi)
return;
clear all
P = setdefaultparams_steak();
return;
phii=0:0.01:1;
for i = 1 : length(phii)
    pi(i) = pi_mix_dim(P,273,phii(i));
    chii(i)= chi(P, 273, phii(i));
end
%plot(phii,chii)


%function y = chi(P, T, phi)
%    y = chi_p( P,T) - (chi_p(P, T)-P.chi_0)*(1-phi)^2;
    %y = .86;
%end



    function x = chi_p(P, T)
        exp_term = P.A*exp(P.gamma*(P.T_0*(ones(size(T))+T*P.nu/P.alpha)-P.T_e*ones(size(T))));
       x =  P.chi_pn*ones(size(T)) + (P.chi_pd-P.chi_pn)*ones(size(T))./(ones(size(T))+exp_term);
    end



function y = pi_mix_dim(P,T,phi)
%    pi_e_dim = P.R.*T*P.rho_s/P.M_c*(phi.^(1/3)*P.phi_0^(2/3)-phi/2);
    pi_e_dim = 4e7/.6/273*T*(phi.^(1/3)*P.phi_0^(2/3)-phi/2);

    pi_m_dim = P.R*T/P.V_f*(log(1-phi)+phi+chi(P,T,phi)*phi^2);
    y = pi_e_dim+pi_m_dim;
end

function phi = root_dim(P,T)
fun = @(x)(pi_mix_dim(P,T,x)); % function
x0 = .6; % initial point
phi = fzero(fun,x0);
end


% Returns root of pi_sw(phi,T) = 0 for given T
% Default tolerance 1e-16
function [x,fval,exitflag,output] = root(P,T,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)(pi_el(P,T,x)+pi_mix(P,T,x)); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end