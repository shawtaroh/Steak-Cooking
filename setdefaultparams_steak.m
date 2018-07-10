%
% 	setdefaultparams_steak.m -- return a data structure with all the parameters
%                       needed to perform a crayfish simulation
%	@project JMU_REU_18
%	@date 6/4/18
%
function P = setdefaultparams_steak
  addpath('./steaksrc/');

    % define the spatial domain
    %
    P.theta = 1.4583; % domain aspect ratio
    P.xmin = 0;
    P.ymin = 0;

    % grid points in each direction
    %

    %  output information
    %
    P.prefix       = 'steak_default';
    P.datadir='.\tmp_3\';
    P.outevery     = 60;  
   
    P.T_0 = 273+7; % Initial Temperature (K)
    P.T_D = 273+225; % Maximum Temperature (K)
    
    % Table 1 Material Properties 
    
    P.rho_f = 1000; % Density of fluid (kg/m^3)
    P.rho_s = 1300; % Density of protein matrix (kg/m^3)
    P.rho_a = 1.225; % Density of air (kg/m^3)

    P.c_f = 4178; % Specific heat of fluid (J / kg.K)
    P.c_s = 2008; % Specific heat of solid (J / kg.K)
    P.c_a = 1006; % Specific heat of air (J / kg.K)
    
    P.k_f = 0.57; % Heat conductivity of fluid (W / m.K)
    P.k_s = 0.18; % Heat conductivity of solid (W / m.K)
    P.k_a = 0.02587; % Heat conductivity of air (W / m.K)
    
    P.r = 2.257e6; % Latent heat of evaporation (J/kg)
    P.h = 15; % Heat transfer coefficient (W/ m^2.K) (Depends on oven, ranges 15-50)
    P.c_0 = 0.24; % Vapor concentration in boundary layer at 100% Humidity at 50 deg C(kg/m^3)
    
    P.D_a = 2.5e-5; % Diffusion coefficient of water vapor in air (m^2/s)
    P.D = P.k_f/P.rho_f/P.c_f; % Thermal diffusivity of fluid (m^2/s)
    P.Diff = 7e-9; % Diffusion coefficient of moisture [Feyesian]
    
    P.kappa_perp = 1e-17; % Permiability of solid in x-direction (m^2)
    P.kappa_par = P.kappa_perp*1.2; % Permiability of solid in y-direction (m^2)
    
    P.M_f = 1.8e-2; % Molar mass of fluid (kg/mol)
    P.V_f = P.M_f/P.rho_f; % Molar volume of fluid (m^3/mol)
    P.M_c = 6; % Elastin molar weight (kg/mol), myosin is 500
    
    
    % Universal Constants
    
    P.R = 8.314; % Gas constant (J/ mol.K)
    
    
    % Fitting Parameters
    
    
    P.chi_0 = 0.5; % Flory-Huggins interaction of fully hydrated polymer
    %P.chi_pn = 0.8; % Flory-Huggins interaction of dry solid protein matrix
    %P.chi_pd = 1.4; % Flory-Huggins interaction of denatured protein matrix
    
    % Alternative Values (Discuss with Shawn)
    %P.chi_pd=.82;
    %P.chi_pn=.675;
    %P.phi_0 = 0.174029030084110;
    % for phi*(T=0) = %17.403, phi*(T=1)=%43.472
    
    
    P.phi_0 = 0.215761926598984;
    P.chi_pd=.9;
    P.chi_pn=.7;
    
    
    %P.chi_pd=.92;
    %P.chi_pn=.72;
    %P.phi_0 = 0.255301981878320;

    
    P.A = 30; % Sigmoidal Fitting
    P.gamma = -0.25; % Sigmoidal Fitting
    P.T_e = 325; % Sigmoidal Fitting
    
    P.X_m = 0.08 ; % GAB a_w moisture content of the first monolayer of absorbed water
    P.p_sat_0 = 597; % Empirical Tetens coefficient (N/m^2)
    
    
        % natural scaling parameters
    %
    P.l = P.k_s/P.h; % length scale (m)
    P.t_0 = P.l^2/P.D; % time scale (s)
    P.mu_0 = mu(P.T_0); % viscosity scale ()
    P.pi_0 = P.mu_0/P.t_0; % pressure scale (N/m^2)
    P.j_evap_0 = P.rho_f*P.l/P.t_0; % evaporative mass flux scale (kg / m^2.s)
    
    
    ly = 0.055/P.l; %ideal thickness
    lx = 0.08/P.l; %ideal width
    P.Ny = 24;
    P.Nx = round(lx/ly*P.Ny);
    P.theta = P.Nx/P.Ny;
    P.xmin = 0;
    P.ymin = 0;
    P.Ly   = sqrt(lx*ly/P.theta); %simulation thickness
    P.Lx   = sqrt(lx*ly*P.theta); %simulation width
    
    % time stepping info
    %
    P.dt = 2e-5;
    P.Nt = round(7200/P.t_0/P.dt);    
    
    
    % Non-dimensional parameters
    
    P.nu = P.rho_s*P.c_s/P.rho_f/P.c_f;
    P.omega = P.k_s/P.k_f;
    P.lambda = P.l^2*P.r*P.rho_f/P.k_s/P.t_0/(P.T_D-P.T_0);
    P.lambda_bar = P.r*P.c_0/(P.T_D-P.T_0)*((P.D_a/P.k_a)^2/P.rho_a/P.c_a)^(1/3);
    P.alpha = P.nu*P.T_0/(P.T_D-P.T_0);
    
    P.N_c = 0.529;
    
    P.beta_el = P.t_0*P.R*P.rho_s*P.T_0/P.mu_0/P.M_c;
    P.beta_el_vm = P.R*P.N_c/P.V_f*P.T_0*P.t_0/P.mu_0;
    P.beta_mix = P.t_0*P.R*P.T_0/P.mu_0/P.V_f;
    
    P.kappa_par_hat = P.kappa_par / P.l^2;
    P.kappa_perp_hat = P.kappa_perp / P.l^2;
    
end