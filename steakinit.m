%
% steak -- return a data structure containing variables needed to
%              solve the NS equations at each time
%
function S=steakinit(P) 


  S.t=0;
  S.w  = zeros(P.Ny,P.Nx,2);        % initial velocity
  %S.phi = root(P,0)*ones(P.Ny,P.Nx); % Initially pi_sw=0
  S.phi = P.phi_0*ones(P.Ny,P.Nx); % Initially pi_sw=0
  % Code for non-uniform initial Phi
    %  S.phi = .17*ones(P.Ny,P.Nx);
    %S.phi(1,:)=P.phi_0*ones(size(S.phi(1,:)));
    %S.phi(:,1)=P.phi_0*ones(size(S.phi(:,1)));
    %S.phi(end,:)=P.phi_0*ones(size(S.phi(end,:)));
    %S.phi(:,end)=P.phi_0*ones(size(S.phi(:,end)));
 
% Enforce Boundary Conditions
  S.T = zeros(P.Ny,P.Nx);
  %h_0 = P.Lx/P.Nx;
  
  % previous time step stuff -- initialize with current values
  %                these will be set after the first time step
  %                
  S.wold    = S.w;
  S.Told = S.T;
  S.phiold    = S.phi;

  % make a flag to indicate a new initialization
  %
  S.isFirst = 1;
  
end
  

% Returns root of pi_sw(phi,T) = 0 for given T
% Default tolerance 1e-16
function phi = root(P,T)
    fun = @(x)(pi_el(P,T,x)+pi_mix(P,T,x)); % function
    x0 = P.phi_0; % initial point
    phi = fzero(fun,x0);
end
