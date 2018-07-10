% nonline_roots_coupled.m - Solves boundary conditions (preprocess)
% for P.T_D=100, P.T_0=0
P = setdefaultparams_steak();
tic
step = 1e-2; % Should use small step
T = 0:step:1;
phi = 0:step:1;

T_left([0 0],P,0,0.4)

for i = 1:length(T)
    for j = 1:length(phi)
        y_left(:,i,j)=T_left([T(i) phi(j)],P,0,.4);
        mod_y(i,j) = y_left(1,i,j)^2+y_left(2,i,j)^2;
    end
end

return;


rootTable = NaN(length(T),length(phi),4,2);
for i = 1:length(T)
    i/length(T)*100
    for j =1:length(phi)
        r_top=root_top(P,T(i),phi(j));
        r_bottom=root_bottom(P,T(i),phi(j));
        r_left=root_left(P,T(i),phi(j));
        r_right=root_right(P,T(i),phi(j));
        rootTable(i,j,:,:) = [r_top; r_bottom; r_left; r_right];
    end
end
toc

save('rootTable.mat','rootTable');

%Returns Left Boundary
function y_left = T_left(x,P,T_int,Phi_int)
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(x(2)).^(-1).*xi(P,x(1)));
    h_int = h_0*sqrt(P.phi_0*(Phi_int).^(-1).*xi(P,T_int));
    T2 = T_int;

    y_left(1) = 1/(P.omega*(1-x(2))+x(2));
    y_left(1) = y_left(1)*(x(1)-T2)/(1/2*H+1/2*h_int)+x(1)-1;
    % j_evap, includes lambda coefficient
    y_left(1) = y_left(1)+j_evap(P,x(1),x(2));
    y_left(2) = (pi_el(P,x(1),x(2))+pi_mix(P,x(1),x(2)));

        % Energy Balance
    % j_evap, includes lambda coefficient
    % pi_sw = 0

end
   
%Returns Left Boundary
function y_right = T_right(x,P,T_int,Phi_int)

    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(x(2)).^(-1).*xi(P,x(2)));
    h_int = h_0*sqrt(P.phi_0*(Phi_int).^(-1).*xi(P,T_int));
    T2 = T_int;

    y_right(1) = 1/(P.omega*(1-x(2))+x(2));
    y_right(1) = y_right(1)*(x(1)-T2)/(1/2*H+1/2*h_int)+x(1)-1;
    % j_evap, includes lambda coefficient
    y_right(1) = y_right(1)+j_evap(P,x(1),x(2));
    y_right(2) = (pi_el(P,x(1),x(2))+pi_mix(P,x(1),x(2)));

end
   
%Returns Top Boundary
% x(1) - T
% x(2) - phi
% Fsolve solves for both s.t. y_top(1)=0 && y_top(2)==0
function y_top = T_top(x,P,T_int,Phi_int)
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(x(2)).^(-1)./xi(P,x(2)));
    h_int = h_0*sqrt(P.phi_0*(Phi_int).^(-1)./xi(P,T_int));
    T2 = T_int;

    y_top(1) = (1-x(2))/P.omega+x(2);
    y_top(1) = y_top(1)*(x(1)-T2)/(1/2*H+1/2*h_int)+x(1)-1;
    % j_evap, includes lambda coefficient
    y_top(1) = y_top(1)+j_evap(P,x(1),x(2));
    y_top(2) = (pi_el(P,x(1),x(2))+pi_mix(P,x(1),x(2)));

end

%Returns Left Boundary
function y_bottom = T_bottom(x,P,T_int,Phi_int)
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(x(2)).^(-1)./xi(P,x(2)));
    h_int = h_0*sqrt(P.phi_0*(Phi_int).^(-1)./xi(P,T_int));
    T2 = T_int;
    
    y_bottom(1) = (1-x(2))/P.omega+x(2);
    y_bottom(1) = y_bottom(1)*(x(1)-T2)/(1/2*H+1/2*h_int)+x(1)-1;
    % j_evap, includes lambda coefficient
    y_bottom(1) = y_bottom(1)+j_evap(P,x(1),x(2));
    y_bottom(2) = (pi_el(P,x(1),x(2))+pi_mix(P,x(1),x(2)));
end
% default fsolve tolerance is 1e-6
function [T phi] = root_bottom(P,T_int,Phi_int)
options = optimset('Display','off','TolFun',1e-6);

fun = @(x)T_bottom(x,P,T_int,Phi_int); % function
x0 = [.2, .5]; % initial point
[T phi] = fsolve(fun,x0,options);
end

function [T phi] = root_top(P,T_int,Phi_int)
options = optimset('Display','off','TolFun',1e-6,'TolX',1e-6);

fun = @(x)T_top(x,P,T_int,Phi_int); % function
x0 = [.2, .5]; % initial point
[T phi] = fsolve(fun,x0,options);
end

function [T phi] = root_left(P,T_int,Phi_int)
options = optimset('Display','off','TolFun',1e-6);

fun = @(x)T_left(x,P,T_int,Phi_int); % function
x0 = [.2, .5]; % initial point
[T phi] = fsolve(fun,x0,options);
end

function [T phi] = root_right(P,T_int,Phi_int)
options = optimset('Display','off','TolFun',1e-6);
fun = @(x)T_right(x,P,T_int,Phi_int); % function
x0 = [.2, .5]; % initial point
[T phi] = fsolve(fun,x0,options);
end
% Default tolerance 1e-16
function phi = root(P,T)
fun = @(x)(pi_el(P,T,x)+pi_mix(P,T,x)); % function
x0 = .6; % initial point
phi = fzero(fun,x0);
end