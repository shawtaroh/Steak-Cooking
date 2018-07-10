%
% steak_sim -- simulates steak
%
%
function steak_sim_altBC( P )
addpath('./steaksrc/');

% output all parameters to a file
%
paramfile( P );

% Make Movie while simulating
writerObj1 = VideoWriter('steak_movie_Evaopration');
writerObj1.FrameRate=5;
open(writerObj1);
Datadir =P.datadir;
% open a figure for plotting
%
scrsz = get(0,'ScreenSize');
fig = figure('Position',[0 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);

% initialize phi s.t. PI_SW = 0
% Optional if not using equilibrium values in setdefaultparams
%P.phi_0 = root(P,0);

outcount = 0;
t=0;
S  = steakinit(P);

h_0 = P.Lx/P.Nx;
h = h_0*ones(P.Ny,P.Nx,2);

% Derivative Operators
% ex. Usage df/dx = (f*cDx)./fH(:,:,1)
fDx = getFDX(P.Nx);
fDy = getFDY(P.Ny);
bDx = getBDX(P.Nx);
bDy = getBDY(P.Ny);

% Main Loop in Time
dt = P.dt;
Nt = P.Nt;
endFlag=0;

output(S,h,t,P.datadir,P.prefix,outcount);
outcount = 1;



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
        
    phi_w_2_fdy = (fDy*((ones(size(S.phi))-S.phiold).*S.w(:,:,2)))./fH(:,:,2);
    phi_w_2_bdy = (bDy*((ones(size(S.phi))-S.phiold).*S.w(:,:,2)))./bH(:,:,2);
    phi_w_2_cdy = 1/2*(phi_w_2_fdy+phi_w_2_bdy);
    % Central Difference divergence
    div_phi_w_cd = phi_w_1_cdx+phi_w_2_cdy;
    % Forward time difference
    %S.phi = max(S.phiold+dt*div_phi_w_cd,0.05*ones(size(S.phi)));
    S.phi = S.phiold+dt*(div_phi_w_cd+P.D_w*P.t_0/P.l^2*(Dxx(S.phiold,h)+Dyy(S.phiold,h)));
    
    
    % Energy Balance
    % 2. Update T
    S.Told = S.T;
    %Finite Differences
    phi_fdx = (S.phiold*fDx)./fH(:,:,1);
    phi_bdx = (S.phiold*bDx)./bH(:,:,1);
    phi_cdx = 1/2*(phi_fdx+phi_bdx);
    
    phi_fdy = (fDy*S.phiold)./fH(:,:,2);
    phi_bdy = (bDy*S.phiold)./bH(:,:,2);
    phi_cdy = 1/2*(phi_fdy+phi_bdy);
    
    T_fdx = (S.T*fDx)./fH(:,:,1);
    T_bdx = (S.T*bDx)./bH(:,:,1);
    T_cdx = 1/2*(T_fdx+T_bdx);
    
    T_fdy = (fDy*S.T)./fH(:,:,2);
    T_bdy = (bDy*S.T)./bH(:,:,2);
    T_cdy = 1/2*(T_fdy+T_bdy);
    
    T_xx = Dxx(S.T,h);
    T_yy = Dyy(S.T,h);
    
    %Isolating T^{k+1} from 7.5.4

    RHS = P.omega^2*ones(size(S.phiold))./(P.omega*(ones(size(S.phiold))-S.phiold)+ones(size(S.phiold))).^2;
    RHS = RHS.*phi_cdx.*T_cdx;
    RHS = RHS+P.omega*ones(size(S.phiold))./(P.omega*(ones(size(S.phiold))-S.phiold)+ones(size(S.phiold))).*T_xx;
    
    RHS = RHS + (P.omega-1)*(phi_cdy.*T_cdy);
    RHS = RHS + (ones(size(S.phiold))-(1-P.omega)*S.phiold).*T_yy;
    
    RHS = RHS + (S.phiold-ones(size(S.phiold))).*(S.w(:,:,1).*T_cdx+S.w(:,:,2).*T_cdy);
    
    RHS = RHS -(P.alpha*ones(size(S.T))+(P.nu*ones(size(S.T))+P.D_w*P.t_0/P.l^2*(Dxx(S.phiold,h)+Dyy(S.phiold,h))).*S.T).*(S.phi-S.phiold)/dt;
    RHS = dt*RHS./(ones(size(S.phiold))-(1-P.nu)*S.phiold);
    S.T = S.T + RHS;
    
    
        % adjust phi on boundary (using old phi and T in j_evap_momentum)
    for j=2:P.Ny-1
        S.phi(j,1)=S.phi(j,2);
        S.phi(j,P.Nx)=S.phi(j,P.Nx-1);
    end
    for i=2:P.Nx-1
        S.phi(1,i)=S.phi(2,i);
        S.phi(P.Ny,i)=S.phi(P.Ny-1,i);
    end
    % take corners to be average of adjacent points
    S.phi(1,1)=(S.phi(1,2)+S.phi(2,1))/2;
    S.phi(P.Ny,1)=(S.phi(P.Ny,2)+S.phi(P.Ny-1,1))/2;
    S.phi(1,P.Nx)=(S.phi(1,P.Nx-1)+S.phi(2,P.Nx))/2;
    S.phi(P.Ny,P.Nx)=(S.phi(P.Ny,P.Nx-1)+S.phi(P.Ny-1,P.Nx))/2;
    
    % ensure phi not to large
    %s.phi=min(S.phi,0.25);

    % set T on boundary (using old T in j_evap)
    for j=2:P.Ny-1
        a=2/(h(j,1,1)+h(j,2,1))/(P.omega*(1-S.phi(j,1))+S.phi(j,1)); %helpful coefficient
        S.T(j,1)=(S.T(j,2)*a+1-j_evap(P,S.Told(j,1),S.phi(j,1),t)) / (1+a);
        a=2/(h(j,P.Nx,1)+h(j,P.Nx-1,1))/(P.omega*(1-S.phi(j,P.Nx))+S.phi(j,P.Nx));
        S.T(j,P.Nx)=(S.T(j,P.Nx-1)*a+1-j_evap(P,S.Told(j,P.Nx),S.phi(j,P.Nx),t)) / (1+a);
    end
    for i=2:P.Nx-1
        a=2/(h(1,i,2)+h(2,i,2))*((1-S.phi(1,i))/P.omega+S.phi(1,i));
        S.T(1,i)=(S.T(2,i)*a+1-j_evap(P,S.Told(1,i),S.phi(1,i),t)) / (1+a);
        a=2/(h(P.Ny,i,2)+h(P.Ny-1,i,2))*((1-S.phi(P.Ny,i))/P.omega+S.phi(P.Ny,i));
        S.T(P.Ny,i)=(S.T(P.Ny-1,i)*a+1-j_evap(P,S.Told(P.Ny,i),S.phi(P.Ny,i),t)) / (1+a);
    end
    % take corners to be average of adjacent points
    S.T(1,1)=(S.T(1,2)+S.T(2,1))/2;
    S.T(P.Ny,1)=(S.T(P.Ny,2)+S.T(P.Ny-1,1))/2;
    S.T(1,P.Nx)=(S.T(1,P.Nx-1)+S.T(2,P.Nx))/2;
    S.T(P.Ny,P.Nx)=(S.T(P.Ny,P.Nx-1)+S.T(P.Ny-1,P.Nx))/2;
    
    %update h on boundaries
%     h(1,:,1)=h_0*sqrt(P.phi_0*(S.phi(1,:)).^(-1).*xi(P,S.T(1,:)));
%     h(1,:,2)=h_0*sqrt(P.phi_0*(S.phi(1,:)).^(-1)./xi(P,S.T(1,:)));
%     h(P.Ny,:,1)=h_0*sqrt(P.phi_0*(S.phi(P.Ny,:)).^(-1).*xi(P,S.T(P.Ny,:)));
%     h(P.Ny,:,2)=h_0*sqrt(P.phi_0*(S.phi(P.Ny,:)).^(-1)./xi(P,S.T(P.Ny,:)));
%     h(:,1,1)=h_0*sqrt(P.phi_0*(S.phi(:,1)).^(-1).*xi(P,S.T(:,1)));
%     h(:,1,2)=h_0*sqrt(P.phi_0*(S.phi(:,1)).^(-1)./xi(P,S.T(:,1)));
%     h(:,P.Nx,1)=h_0*sqrt(P.phi_0*(S.phi(:,P.Nx)).^(-1).*xi(P,S.T(:,P.Nx)));
%     h(:,P.Nx,2)=h_0*sqrt(P.phi_0*(S.phi(:,P.Nx)).^(-1)./xi(P,S.T(:,P.Nx)));


    %   Shrinkage
    % 3. Update h on bulk
    
%     h(2:end-1,2:end-1,2) = h_0*sqrt(P.phi_0*(S.phi(2:end-1,2:end-1)).^(-1)./xi(P,S.T(2:end-1,2:end-1)));
%     h(2:end-1,2:end-1,1) = h_0*sqrt(P.phi_0*(S.phi(2:end-1,2:end-1)).^(-1).*xi(P,S.T(2:end-1,2:end-1)));
    phi_avg=mean(mean(S.phi));
    xi_avg=mean(mean(xi(P,S.T)));
    h(:,:,1)=h_0*(P.g*sqrt(P.phi_0/phi_avg*xi_avg)*ones(P.Ny,P.Nx)+(1-P.g)*sqrt(P.phi_0*(S.phi).^(-1).*xi(P,S.T)));
    h(:,:,2)=h_0*(P.g*sqrt(P.phi_0/phi_avg/xi_avg)*ones(P.Ny,P.Nx)+(1-P.g)*sqrt(P.phi_0*(S.phi).^(-1)./xi(P,S.T)));


    % Update finite difference H's
    fH = getFH(P.Ny,P.Nx,h);
    bH = getBH(P.Ny,P.Nx,h);
    
    % Momentum Balance (7.5.2)
    % Pseudo-code Step 6
        
    S.wold = S.w;
    
    PI_sw = pi_el(P,S.T,S.phi) + pi_mix(P,S.T,S.phi);
    PI_sw(1,:) = zeros(size(PI_sw(1,:)));
    PI_sw(:,1) = zeros(size(PI_sw(:,1)));
    PI_sw(end,:) = zeros(size(PI_sw(end,:)));
    PI_sw(:,end) = zeros(size(PI_sw(:,end)));
    % Finite Difference Derivatives
    pi_fdx = (PI_sw*fDx)./fH(:,:,1);
    pi_bdx = (PI_sw*bDx)./bH(:,:,1);
    pi_cdx = 1/2*(pi_fdx+pi_bdx);
    
    pi_fdy = (fDy*PI_sw)./fH(:,:,2);
    pi_bdy = (bDy*PI_sw)./bH(:,:,2);
    pi_cdy = 1/2*(pi_fdy+pi_bdy);
    
   % Momentum Balance (7.5.2)
    S.w(:,:,1)=-P.kappa_perp_hat*pi_cdx./mu_hat(P,S.T)./(ones(size(S.phi))-S.phi);
    S.w(:,:,2)=-P.kappa_par_hat*pi_cdy./mu_hat(P,S.T)./(ones(size(S.phi))-S.phi);
    
    
    % Stop for internal temperature 65 deg C
    if(S.T(floor(P.Ny/2),floor(P.Nx/2))*(P.T_D-P.T_0)+P.T_0-273>65)
        %endFlag = 1;
        %disp("Dinner's ready");
    end
    
    S.isFirst = 0;
    %% Video and Saving Files
    
    if( mod(k,P.outevery)==0 )
        %if( false )
        
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
        %hpc = pcolor(x,y,S.phi);
        set(hpc,'edgecolor','none','facelighting','flat','facecolor','interp');
        
        caxis([0 (65+273-P.T_0)/(P.T_D-P.T_0)])
        %caxis([0.17 0.5])
        
        grid on;
        colorbar;
        scale = 1.25;
        axis([-scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2 -scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2]);
        % Plot Grid of x,y
        for i=1:Nx-1
            for j=1:Ny-1
                plot([x(j,i) x(j,i+1)], [y(j,i) y(j,i+1)],'r', 'LineWidth',.5);% right line
                plot([x(j,i) x(j+1,i)], [y(j,i) y(j+1,i)], 'r','LineWidth',.5);% top line
                if(i==Nx-1)
                    plot([x(j,i+1) x(j+1,i+1)],[y(j,i+1) y(j+1,i+1)],'r', 'LineWidth',.5);
                end
            end
            plot([x(Ny,i) x(Ny,i+1)],[y(Ny,i) y(Ny,i+1)],'r', 'LineWidth',.5);
        end
        %Plot Velocity Vector Field
        
        quiver(x,y,S.w(:,:,1),S.w(:,:,2),'b','LineWidth',.5);
        ylabel('(cm)');
        xlabel(sprintf('Normalized time = %.4f, Real Time = %.4f (sec)',t,t*P.t_0),'FontSize',12)
        
        
        hold off
        frame=getframe(fig);
        
        writeVideo(writerObj1,frame);
        %end
        
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
        %fprintf('Saving file %s \n',filename);
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
    fDx = diag(ones(Nx-1,1),-1)-diag(ones(Nx,1));
    fDx(Nx-1,Nx)=-1;
    fDx(Nx,Nx)=1;
end

% returns forward 1st y-derivative operator (7.5)
function fDy = getFDY(Ny)
    fDy = diag(ones(Ny-1,1),1)-diag(ones(Ny,1));
    fDy(Ny,Ny-1)=-1;
    fDy(Ny,Ny)=1;
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
    bDx = diag(ones(Nx,1))-diag(ones(Nx-1,1),1);
    bDx(1,1)=-1;
    bDx(2,1)=1;
end

% returns backward 1st y-derivative operator (7.5)
function bDy = getBDY(Ny)
    bDy = diag(ones(Ny,1))-diag(ones(Ny-1,1),-1);
    bDy(1,1)=-1;
    bDy(1,2)=1;
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