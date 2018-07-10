%
% steak_sim -- simulates steak
%
%
function steak_sim_fastest_VDS( P )
addpath('./steaksrc/');
load('low_phi_root_interp_175.mat','T_interp','phi_interp')
% output all parameters to a file
%
paramfile( P );

% Make Movie while simulating
writerObj1 = VideoWriter('steak_movie_altBC');
writerObj1.FrameRate=5;
open(writerObj1);
Datadir =P.datadir;
% open a figure for plotting
%
scrsz = get(0,'ScreenSize')
fig = figure('Position',[scrsz(3)/2 scrsz(4)/2 720 720])

% initialize phi s.t. PI_SW = 0
% Optional if not using equilibrium values in setdefaultparams
%P.phi_0 = root(P,0);

outcount = 0;
t=0;
S  = steakinit(P);

h_0 = P.Lx/P.Nx;
h = h_0*ones(P.Ny,P.Nx,2);

% Main Loop in Time
dt = P.dt;
Nt = P.Nt;
endFlag=0;

output(S,h,t,P.datadir,P.prefix,outcount);
outcount = 1;

% Derivative Operators
% ex. Usage df/dx = (f*cDx)./fH(:,:,1)
fDx = getFDX(P.Nx);
fDy = getFDY(P.Ny);
bDx = getBDX(P.Nx);
bDy = getBDY(P.Ny);

fH = getFH(P.Ny,P.Nx,h);
bH = getBH(P.Ny,P.Nx,h);

for k=1:Nt
    
    tic;
    % compute the time
    %
    t = k*dt;
    S.t=t;
    %   Continuity
    % 1.  Phi^{k+1} <- Phi^{k}, T^{k}, w^{k}, h^{k}
    
    S.phiold=S.phi;
    % Finite Differences
    phi_w_1_fdx = (((ones(size(S.phi))-S.phiold).*S.w(:,:,1))*fDx)./fH(:,:,1);
    phi_w_1_bdx = (((ones(size(S.phi))-S.phiold).*S.w(:,:,1))*bDx)./bH(:,:,1);
    phi_w_1_cdx = 1/2*(phi_w_1_fdx+phi_w_1_bdx);
    phi_w_1_cdx(1:P.Ny,1) = phi_w_1_fdx(1:P.Ny,1);
    phi_w_1_cdx(1:P.Ny,P.Nx) = phi_w_1_bdx(1:P.Ny,P.Nx);
        
    phi_w_2_fdy = (fDy*((ones(size(S.phi))-S.phiold).*S.w(:,:,2)))./fH(:,:,2);
    phi_w_2_bdy = (bDy*((ones(size(S.phi))-S.phiold).*S.w(:,:,2)))./bH(:,:,2);
    phi_w_2_cdy = 1/2*(phi_w_2_fdy+phi_w_2_bdy);
    phi_w_2_cdy(1, 1:P.Nx) = phi_w_2_fdy(1, 1:P.Nx);
    phi_w_2_cdy(P.Ny, 1:P.Nx) = phi_w_2_bdy(P.Ny, 1:P.Nx);
    % Central Difference divergence
    div_phi_w_cd = phi_w_1_cdx+phi_w_2_cdy;
    % Diffusion
    phi_xx = Dxx(S.phi,h);
    phi_yy = Dyy(S.phi,h);
    diffusion = P.t_0*P.Diff/P.l^2*(phi_xx+phi_yy);
    % Forward time difference
    S.phi = S.phiold+dt*(div_phi_w_cd+diffusion);

    
    % Energy Balance
    % 2. Update T
    S.Told = S.T;
    %Finite Differences
    phi_fdx = (S.phiold*fDx)./fH(:,:,1);
    phi_bdx = (S.phiold*bDx)./bH(:,:,1);
    phi_cdx = 1/2*(phi_fdx+phi_bdx);
    phi_cdx(1:P.Ny,1) = phi_fdx(1:P.Ny,1);
    phi_cdx(1:P.Ny,P.Nx) = phi_bdx(1:P.Ny,P.Nx);
    
    phi_fdy = (fDy*S.phiold)./fH(:,:,2);
    phi_bdy = (bDy*S.phiold)./bH(:,:,2);
    phi_cdy = 1/2*(phi_fdy+phi_bdy);
    phi_cdy(1, 1:P.Nx) = phi_fdy(1, 1:P.Nx);
    phi_cdy(P.Ny, 1:P.Nx) = phi_bdy(P.Ny, 1:P.Nx);
    
    T_fdx = (S.T*fDx)./fH(:,:,1);
    T_bdx = (S.T*bDx)./bH(:,:,1);
    T_cdx = 1/2*(T_fdx+T_bdx);
    T_cdx(1:P.Ny,1) = T_fdx(1:P.Ny,1);
    T_cdx(1:P.Ny,P.Nx) = T_bdx(1:P.Ny,P.Nx);
    
    T_fdy = (fDy*S.T)./fH(:,:,2);
    T_bdy = (bDy*S.T)./bH(:,:,2);
    T_cdy = 1/2*(T_fdy+T_bdy);
    T_cdy(1, 1:P.Nx) = T_fdy(1, 1:P.Nx);
    T_cdy(P.Ny, 1:P.Nx) = T_bdy(P.Ny, 1:P.Nx);
    
    T_xx = Dxx(S.T,h);
    T_yy = Dyy(S.T,h);
    
    %Isolating T^{k+1} from 7.5.4

    RHS = P.omega^2*ones(size(S.phiold))./(P.omega*(ones(size(S.phiold))-S.phiold)+ones(size(S.phiold))).^2;
    RHS = RHS.*phi_cdx.*T_cdx;
    RHS = RHS+P.omega*ones(size(S.phiold))./(P.omega*(ones(size(S.phiold))-S.phiold)+ones(size(S.phiold))).*T_xx;
    
    RHS = RHS + (P.omega-1)*(phi_cdy.*T_cdy);
    RHS = RHS + (ones(size(S.phiold))-(1-P.omega)*S.phiold).*T_yy;
    
    RHS = RHS + (S.phiold-ones(size(S.phiold))).*(S.w(:,:,1).*T_cdx+S.w(:,:,2).*T_cdy);
    
    RHS = RHS -(P.alpha*ones(size(S.T))+P.nu*S.T).*(S.phi-S.phiold)/dt;
    RHS = RHS + P.Diff/P.D*(phi_xx+phi_yy).*(S.T+P.alpha/P.nu*ones(size(S.T)));
    RHS = dt*RHS./(ones(size(S.phiold))-(1-P.nu)*S.phiold);
    S.T = S.T + RHS;
    
    %   Shrinkage
    % 3. Update h on bulk
    
    h(2:end-1,2:end-1,2) = h_0*(P.phi_0*(S.phi(2:end-1,2:end-1)).^(-1)./xi(P,S.T(2:end-1,2:end-1))).^(1/3);
    h(2:end-1,2:end-1,1) = h_0*(P.phi_0*(S.phi(2:end-1,2:end-1)).^(-1).*xi(P,S.T(2:end-1,2:end-1))).^(1/3);
    
    % Average the corner temperatures (Corner Temperatures don't actually do anything, but the phi at the corners matter b/c of momentum equation gradient pi_sw)
    % These corners might not matter at all

    % Update finite difference H's
    fH = getFH(P.Ny,P.Nx,h);
    bH = getBH(P.Ny,P.Nx,h);
    
    % Momentum Balance (7.5.2)
    % Pseudo-code Step 6
        
    S.wold = S.w;
    S.T(S.T<0)=0;
    phi_eq = interp1(T_interp,phi_interp,min(ones(size(S.T)),S.T*(P.T_D-273)/175)); % Use previous T on boundary to interpolate preprocessed root data
    n_eq = (ones(size(phi_eq))-phi_eq)./(ones(size(phi_eq))+.3*phi_eq);
    n = (ones(size(S.phi))-S.phi)./(ones(size(S.phi))+.3*S.phi);
    
    
    PI_sw = (n-n_eq);

    PI_sw(1,:) = zeros(size(PI_sw(1,:)));
    PI_sw(:,1) = zeros(size(PI_sw(:,1)));
    PI_sw(end,:) = zeros(size(PI_sw(end,:)));
    PI_sw(:,end) = zeros(size(PI_sw(:,end)));
    % Finite Difference Derivatives
    pi_fdx = (PI_sw*fDx)./fH(:,:,1);
    pi_bdx = (PI_sw*bDx)./bH(:,:,1);
    pi_cdx = 1/2*(pi_fdx+pi_bdx);
    pi_cdx(1:P.Ny,1) = pi_fdx(1:P.Ny,1);
    pi_cdx(1:P.Ny,P.Nx) = pi_bdx(1:P.Ny,P.Nx);
    
    pi_fdy = (fDy*PI_sw)./fH(:,:,2);
    pi_bdy = (bDy*PI_sw)./bH(:,:,2);
    pi_cdy = 1/2*(pi_fdy+pi_bdy);
    pi_cdy(1, 1:P.Nx) = pi_fdy(1,1:P.Nx);
    pi_cdy(P.Ny, 1:P.Nx) = pi_bdy(P.Ny,1:P.Nx);
    
   % Momentum Balance (7.5.2)
    S.w(:,:,1)=-10^-7*P.t_0/P.l^2*pi_cdx./(ones(size(S.phi))-S.phi);
    S.w(:,:,2)=-10^-7*P.t_0/P.l^2*pi_cdy./(ones(size(S.phi))-S.phi);
    
    
    % Coupled Nonlinear Boundary Conditions
    % 4. Update Boundary Conditions
    % 5. Shrinkage on boundary
    %Top/Bottom T & Phi
    for i=2:P.Nx-1
            S.phi(1,i)= interp1(T_interp,phi_interp,min(1,S.Told(1,i)*(P.T_D-273)/175)); % Use previous T on boundary to interpolate preprocessed root data
            S.phi(P.Ny,i)= interp1(T_interp,phi_interp,min(1,S.Told(P.Ny,i)*(P.T_D-273)/175)); % Use previous T on boundary to interpolate preprocessed root data
            %S.phi(1,i) = S.phi(2,i);
            %S.phi(P.Ny,i) = S.phi(P.Ny-1,i);
            if(S.isFirst)
                [r_top,fval,exitflag,out]=root_top(P,S,h,i,.1);
                if(exitflag>0 && imag(r_top(1)) ==0 && imag(fval(1)) ==0)
                    S.T(1,i)=r_top(1);
                else
                    disp('ROOT SOLVING ERROR \n'); 
                    return;
                end
                [r_bot,fval,exitflag,out] = root_bottom(P,S,h,i,.1);
                 if(exitflag>0 && imag(r_bot(1)) ==0 && imag(fval(1)) ==0)
                    S.T(P.Ny,i) = r_bot(1);
                 else
                    disp('ROOT SOLVING ERROR \n'); 
                    return;
                 end
            else
                	a=2/(h(1,i,2)+h(2,i,2))*((1-S.phi(1,i))/P.omega+S.phi(1,i));
                    S.T(1,i)=max(S.Told(1,i),(S.T(2,i)*a+1-j_evap(P,S.Told(1,i),S.phi(1,i),t)) / (1+a));
                    a=2/(h(P.Ny,i,2)+h(P.Ny-1,i,2))*((1-S.phi(P.Ny,i))/P.omega+S.phi(P.Ny,i));
                    S.T(P.Ny,i)=max(S.Told(P.Ny,i),(S.T(P.Ny-1,i)*a+1-j_evap(P,S.Told(P.Ny,i),S.phi(P.Ny,i),t)) / (1+a));
            end
        h(1,i,1) = h_0*sqrt(P.phi_0*(S.phi(1,i)).^(-1).*xi(P,S.T(1,i)));
        h(1,i,2) = h_0*sqrt(P.phi_0*(S.phi(1,i)).^(-1)./xi(P,S.T(1,i)));
        h(P.Ny,i,1) = h_0*sqrt(P.phi_0*(S.phi(P.Ny,i)).^(-1).*xi(P,S.T(P.Ny,i)));
        h(P.Ny,i,2) = h_0*sqrt(P.phi_0*(S.phi(P.Ny,i)).^(-1)./xi(P,S.T(P.Ny,i)));
    end
    %Left/Right T & Phi
    for j=2:P.Ny-1
        S.phi(j,1)= interp1(T_interp,phi_interp,min(1,S.Told(j,1)*(P.T_D-273)/175)); % Use previous T on boundary to interpolate preprocessed root data
        S.phi(j,P.Nx)= interp1(T_interp,phi_interp,min(1,S.Told(j,P.Nx)*(P.T_D-273)/175)); % Use previous T on boundary to interpolate preprocessed root data
        %    S.phi(j,1) = S.phi(j,2);
        %    S.phi(j,P.Nx) = S.phi(j,P.Nx-1);
        if(S.isFirst)
            [r_left,fval,exitflag,out]=root_left(P,S,h,j,.1);
            if(exitflag>0 && imag(r_left(1)) ==0 && imag(fval(1)) ==0 )
                S.T(j,1) = r_left(1);
            else
                disp('ROOT FINDING ERROR\n'); 
                return;
            end
            [r_right,fval,exitflag,out] = root_right(P,S,h,j,.1);
            if(exitflag>0 && imag(r_right(1)) ==0 && imag(fval(1)) ==0)
                S.T(j,P.Nx) = r_right(1);
            else
                disp('ROOT FINDING ERROR\n'); 
                return;
            end
        else
            % set T on boundary (using old T in j_evap)
            a=2/(h(j,1,1)+h(j,2,1))/(P.omega*(1-S.phi(j,1))+S.phi(j,1)); %helpful coefficient
            S.T(j,1)=max(S.Told(j,1),(S.T(j,2)*a+1-j_evap(P,S.Told(j,1),S.phi(j,1),t)) / (1+a));
            a=2/(h(j,P.Nx,1)+h(j,P.Nx-1,1))/(P.omega*(1-S.phi(j,P.Nx))+S.phi(j,P.Nx));
            S.T(j,P.Nx)=max(S.Told(j,P.Nx),(S.T(j,P.Nx-1)*a+1-j_evap(P,S.Told(j,P.Nx),S.phi(j,P.Nx),t)) / (1+a));
 
        end
        h(j,1,1) = h_0*sqrt(P.phi_0*(S.phi(j,1)).^(-1).*xi(P,S.T(j,1)));
        h(j,1,2) = h_0*sqrt(P.phi_0*(S.phi(j,1)).^(-1)./xi(P,S.T(j,1)));
        h(j,P.Nx,1) = h_0*sqrt(P.phi_0*(S.phi(j,P.Nx)).^(-1).*xi(P,S.T(j,P.Nx)));
        h(j,P.Nx,2) = h_0*sqrt(P.phi_0*(S.phi(j,P.Nx)).^(-1)./xi(P,S.T(j,P.Nx)));
    end
   
S.phi(1,1)=(S.phi(1,2)+S.phi(2,1))/2;
S.phi(P.Ny,1)=(S.phi(P.Ny,2)+S.phi(P.Ny-1,1))/2;
S.phi(1,P.Nx)=(S.phi(1,P.Nx-1)+S.phi(2,P.Nx))/2;
S.phi(P.Ny,P.Nx)=(S.phi(P.Ny,P.Nx-1)+S.phi(P.Ny-1,P.Nx))/2;

% take corners to be average of adjacent points
    S.T(1,1)=(S.T(1,2)+S.T(2,1))/2;
    S.T(P.Ny,1)=(S.T(P.Ny,2)+S.T(P.Ny-1,1))/2;
    S.T(1,P.Nx)=(S.T(1,P.Nx-1)+S.T(2,P.Nx))/2;
    S.T(P.Ny,P.Nx)=(S.T(P.Ny,P.Nx-1)+S.T(P.Ny-1,P.Nx))/2;

    % Update h for all sides
    
    %Bottom
    
    h(P.Ny,:,1) = h_0*(P.phi_0*(min(S.phi(P.Ny,:),0.4347*ones(size(S.phi(P.Ny,:))))).^(-1).*xi(P,S.T(P.Ny,:))).^(1/3);
    h(P.Ny,:,2) = h_0*(P.phi_0*(min(S.phi(P.Ny,:),0.4347*ones(size(S.phi(P.Ny,:))))).^(-1)./xi(P,S.T(P.Ny,:))).^(1/3);
    
    %Top
    h(1,:,1) = h_0*(P.phi_0*(min(S.phi(1,:),0.4347*ones(size(S.phi(1,:))))).^(-1).*xi(P,S.T(1,:))).^(1/3);
    h(1,:,2) = h_0*(P.phi_0*(min(S.phi(1,:),0.4347*ones(size(S.phi(1,:))))).^(-1)./xi(P,S.T(1,:))).^(1/3);
    
    %Left
    h(:,1,1) = h_0*(P.phi_0*(min(S.phi(:,1),0.4347*ones(size(S.phi(:,1))))).^(-1).*xi(P,S.T(:,1))).^(1/3);
    h(:,1,2) = h_0*(P.phi_0*(min(S.phi(:,1),0.4347*ones(size(S.phi(:,1))))).^(-1)./xi(P,S.T(:,1))).^(1/3);
    %Right
    h(:,P.Nx,1) = h_0*(P.phi_0*(min(S.phi(:,P.Nx),0.4347*ones(size(S.phi(:,P.Nx))))).^(-1).*xi(P,S.T(:,P.Nx))).^(1/3);
    h(:,P.Nx,2) = h_0*(P.phi_0*(min(S.phi(:,P.Nx),0.4347*ones(size(S.phi(:,P.Nx))))).^(-1)./xi(P,S.T(:,P.Nx))).^(1/3);

    
    % Stop for internal temperature 65 deg C
    if(S.T(floor(P.Ny/2),floor(P.Nx/2))*(P.T_D-P.T_0)+P.T_0-273>65)
        %endFlag = 1;
    end
    
    S.isFirst = 0;
    %% Video and Saving Files
    
    if( mod(k,P.outevery)==0 )
        
        if( mod(k,P.outevery*100)==0 )
        plot(1);
        % Get x,y coordinates from h
        [x y] = visualize(h);
         x = P.l*100*x; % Nondimensional scaling in cm
         y = P.l*100*y;
        hold on;
        Nx = P.Nx;
        Ny = P.Ny;
        % Color Map
        hpc = pcolor(x,y,S.T);
        %    hpc = pcolor(x,y,S.phi);
        set(hpc,'edgecolor','none','facelighting','flat','facecolor','interp');
        
        caxis([0 (65+273-P.T_0)/(P.T_D-P.T_0)])
        
        grid on;
        colorbar;
        scale = 1.25
        axis([-scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2 -scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2]);
        % Plot Grid of x,y
       for i=1:Nx-1
                plot(x(:,i), y(:,i),'r', 'LineWidth',.5);% right line
            for j=1:Ny-1
                if(i==Nx-1)
                    plot([x(j,i+1) x(j+1,i+1)],[y(j,i+1) y(j+1,i+1)],'r', 'LineWidth',.5);
                end
            end
            plot([x(Ny,i) x(Ny,i+1)],[y(Ny,i) y(Ny,i+1)],'r', 'LineWidth',.5);
        end
        
        for j=1:Ny-1
              plot(x(j,:), y(j,:),'r', 'LineWidth',.5);% top line 
        end
        %Plot Velocity Vector Field
        
        quiver(x,y,S.w(:,:,1),S.w(:,:,2),'b','LineWidth',.5);
        ylabel('(cm)');
        xlabel(sprintf('Normalized time = %.4f, Real Time = %.4f (sec)',t,t*P.t_0),'FontSize',12)
        
        
        hold off
        frame=getframe(fig);
        
        writeVideo(writerObj1,frame);
        end
        k
        
        output(S,h,t,P.datadir,P.prefix,outcount);
        outcount = outcount + 1;
        if(endFlag==1)
            error('Breaking out of function');
        end
    end
end

disp("ending")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output data to a file
%
    function output(S,h,t,datadir,prefix,outcount)
        filename = sprintf('%s%s.%5.5i.mat',datadir,prefix,outcount);
        save(filename,'h','S','t');
        fprintf('Saving file %s \n',filename);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output all of the parameters
%
function paramfile( P )
        
        % write an ascii file
        %
        filename = sprintf('%s/%s.param',P.datadir,P.prefix);
        fid = fopen(filename,'wt');
        pfnames = fieldnames(P);
        pfnameschar = char( pfnames );
        for k=1:length(pfnames)
            disp(k)
            field = getfield(P,pfnames{k});
            fprintf(fid,'%s ',pfnameschar(k,:));
            
            if( ischar( field ) )
                fprintf(fid,'     %s\n',field);
            else
                for m=1:length(field)
                    fprintf(fid,'      %g\n',field(m));
                end
            end
        end
        fclose(fid);
        
        % write a mat file for easier processing
        %
        filename = sprintf('%s%s.param.mat',P.datadir,P.prefix);
        save(filename,'P');
    end

    close(writerObj1);

end

% returns forward 1st x-derivative operator (7.5)
function fDx = getFDX(Nx)
    fDx = zeros(Nx);
    for i=1:Nx-1
        fDx(i,i)=-1;
        fDx(i+1,i)=1;
    end
    fDx(Nx,Nx)=-1;
end

% returns forward 1st y-derivative operator (7.5)
function fDy = getFDY(Ny)
    fDy = zeros(Ny);
    for j=1:Ny-1
        fDy(j,j)=-1;
        fDy(j,j+1)=1;
    end
    fDy(Ny,Ny)=-1;
end

% returns denominator H for use with fDx and fDy (7.5)
function fH = getFH(Ny,Nx,h)
    fH = zeros(Ny,Nx,2);
    % x-direction
    for i=1:Nx-1
        for j=1:Ny
           fH(j,i,1) = 0.5*(h(j,i,1)+h(j,i+1,1));
        end
    end
    for j=1:Ny
        fH(j,Nx,1)=h(j,Nx,1);
    end
    % y-direction
    for j=1:Ny-1
        for i=1:Nx
            fH(j,i,2)=0.5*(h(j,i,2)+h(j+1,i,2));
        end
    end
    for i=1:Nx
        fH(Ny,i,2)=h(Ny,i,2);
    end
end

% returns backward 1st x-derivative operator (7.5)
function bDx = getBDX(Nx)
    bDx = zeros(Nx);
    for i=2:Nx
        bDx(i-1,i)=-1;
        bDx(i,i)=1;
    end
    bDx(1,1)=1;
end

% returns backward 1st y-derivative operator (7.5)
function bDy = getBDY(Ny)
    bDy = zeros(Ny);
    for j=2:Ny
        bDy(j,j-1)=-1;
        bDy(j,j)=1;
    end
    bDy(1,1)=1;
end

% returns denominator H for use with bDx and bDy (7.5)
function bH = getBH(Ny,Nx,h)
    bH = zeros(Ny,Nx,2);
    % x-direction
    for i=2:Nx
        for j=1:Ny
            bH(j,i,1) = 0.5*(h(j,i-1,1)+h(j,i,1));
        end
    end
    for j=1:Ny
        bH(j,1,1)=h(j,1,1);
    end
% y-direction
    for j=2:Ny
        for i=1:Nx
            bH(j,i,2)=0.5*(h(j-1,i,2)+h(j,i,2));
        end
    end
    for i=1:Nx
        bH(1,i,2)=h(1,i,2);
    end
end

% Returns Discrete Second Derivative of f given spacing h (7.5)
function D = Dxx(f,h)
    [Ny Nx] = size(f);
    D=zeros(size(f));
    for j=1:Ny
        for i=2:Nx-1
            D(j,i)=2*((f(j,i+1)-f(j,i))/(h(j,i+1,1)+h(j,i,1))-(f(j,i)-f(j,i-1))/(h(j,i,1)+h(j,i-1,1)))/h(j,i,1);
        end
        D(j,1)=D(j,2);
        D(j,Nx)=D(j,Nx-1);
    end
end

% Returns Discrete Second Derivative of f given spacing h (7.5)
function D = Dyy(f,h)
    [Ny Nx] = size(f);
    D=zeros(size(f));
    for i=1:Nx
        for j=2:Ny-1
            D(j,i)=2*((f(j+1,i)-f(j,i))/(h(j+1,i,2)+h(j,i,2))-(f(j,i)-f(j-1,i))/(h(j,i,2)+h(j-1,i,2)))/h(j,i,2);
        end
        D(1,i)=D(2,i);
        D(Ny,i)=D(Ny-1,i);
    end
end

%Returns Nonlinear System for Left Boundary (7.6)
% x(1) - T, x(2) - phi
% Fsolve solves for both s.t. y_left(1)=0 && y_left(2)==0
function y_left = T_left(x,P,S,h,j)
    % Energy Balance
    y_left(1) = 1/(P.omega*(1-S.phi(j,1))+S.phi(j,1));
    T2 = S.T(j,2);
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(S.phi(j,1)).^(-1).*xi(P,x(1)));
    y_left(1) = y_left(1)*(x(1)-T2)/(1/2*H+1/2*h(j,2,1))+x(1)-1;
    % j_evap, includes lambda coefficient
    y_left(1) = y_left(1)+j_evap(P,x(1),S.phi(j,1),S.t);
    % pi_sw = 0
end

%Returns Nonlinear System for Right Boundary (7.6)
% x(1) - T, x(2) - phi
% Fsolve solves for both s.t. y_right(1)=0 && y_right(2)==0
function y_right = T_right(x,P,S,h,j)
    % Energy Balance
    y_right(1) = 1/(P.omega*(1-S.phi(j,end))+S.phi(j,end));
    T2 = S.T(j,end-1);
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(S.phi(j,end)).^(-1).*xi(P,x(1)));
    y_right(1) = y_right(1)*(x(1)-T2)/(1/2*H+1/2*h(j,end-1,1))+x(1)-1;
    % j_evap, includes lambda coefficient
    y_right(1) = y_right(1)+j_evap(P,x(1),S.phi(j,end),S.t);
    % pi_sw = 0
end

%Returns Nonlinear System for Top Boundary (7.6)
% x(1) - T
% x(2) - phi
% Fsolve solves for both s.t. y_top(1)=0 && y_top(2)==0
function y_top = T_top(x,P,S,h,i)
    % Energy Balance
    y_top(1) = (1-S.phi(1,i))/P.omega+S.phi(1,i);
    T2 = S.T(2,i);
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(S.phi(1,i)).^(-1)./xi(P,x(1)));
    y_top(1) = y_top(1)*(x(1)-T2)/(1/2*H+1/2*h(2,i,2))+x(1)-1;
    %j_evap includes lambda bar
    y_top(1) = y_top(1)+j_evap(P,x(1),S.phi(1,i),S.t);
    % pi_sw = 0
end

%Returns Nonlinear System for Bottom Boundary (7.6)
% x(1) - T, x(2) - phi
% Fsolve solves for both s.t. y_bottom(1)=0 && y_bottom(2)==0
function y_bottom = T_bottom(x,P,S,h,i)
    % Energy Balance
    y_bottom(1) = (1-S.phi(end,i))/P.omega+S.phi(end,i);
    T2 = S.T(end-1,i);
    h_0 = P.Lx/P.Nx;
    H = h_0*sqrt(P.phi_0*(S.phi(end,i)).^(-1)./xi(P,x(1)));
    y_bottom(1) = y_bottom(1)*(x(1)-T2)/(1/2*H+1/2*h(end-1,i,2))+x(1)-1;
    % j_evap, includes lambda coefficient
    y_bottom(1) = y_bottom(1)+j_evap(P,x(1),S.phi(end,i),S.t);
    % pi_sw = 0
end

%Returns Roots of Nonlinear System for Bottom Boundary (7.6)
% Fsolve solves for both s.t. y_bottom(1)=0 && y_bottom(2)==0
% default fsolve tolerance is 1e-6
function [x,fval,exitflag,output] = root_bottom(P,S,h,i,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)T_bottom(x,P,S,h,i); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end

%Returns Roots of Nonlinear System for Top Boundary (7.6)
% Fsolve solves for both s.t. y_top(1)=0 && y_top(2)==0
% default fsolve tolerance is 1e-6
function [x,fval,exitflag,output] = root_top(P,S,h,i,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)T_top(x,P,S,h,i); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end

%Returns Roots of Nonlinear System for Left Boundary (7.6)
% Fsolve solves for both s.t. y_left(1)=0 && y_left(2)==0
% default fsolve tolerance is 1e-6
function [x,fval,exitflag,output] = root_left(P,S,h,j,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)T_left(x,P,S,h,j); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end

%Returns Roots of Nonlinear System for Right Boundary (7.6)
% Fsolve solves for both s.t. y_right(1)=0 && y_right(2)==0
% default fsolve tolerance is 1e-6
function [x,fval,exitflag,output] = root_right(P,S,h,j,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)T_right(x,P,S,h,j); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end

% Returns root of pi_sw(phi,T) = 0 for given T
% Default tolerance 1e-16
function [x,fval,exitflag,output] = root(P,S,j,i,x0)
    options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    fun = @(x)(pi_el(P,S.Told(j,i),x)+pi_mix(P,S.Told(j,i),x)); % function
    [x,fval,exitflag,output] = fzero(fun,x0,options);
end