%
% make a movie 
%
Datadir ='.\perm1e-16_nodiff\'
% open a figure for plotting
%
scrsz = get(0,'ScreenSize')
fig = figure('Position',[scrsz(3)/8 scrsz(4)/2 640 640])

% video object
%
infilebase =  sprintf('%s/steak_default',Datadir);

    % load grid and parameters
    %
    paramfile = sprintf('%s.param.mat',infilebase);
    load(paramfile);
    
writerObj1 = VideoWriter('perm1e-16_nodiff_p');
writerObj1.FrameRate=5;
open(writerObj1);

filePattern = fullfile(Datadir, 'steak_default.0*.mat');
matFiles = dir(filePattern);
names = matFiles.name;
%length(matFiles)

% w=[0 0];
% for k_t=1:length(matFiles)/2
%     disp(k_t)
%         % load the grid and parameters
%     load(sprintf('%s%s',Datadir,matFiles(k_t).name));
%     for i =1:P.Nx
%         for j = 1:P.Ny
%             if(S.w(j,i,1)^2+S.w(j,i,2)^2>w(1)^2+w(2)^2)
%                w(1)= S.w(j,i,1);
%                w(2)= S.w(j,i,2);
%                
%             end
%         end
%     end
% end

% (1-phi)w

phi_w=[0 0];
T_max = 0;
phi_max = 0;
phi_min = 1;
for k_t=1:100:length(matFiles)
    disp(k_t)
        % load the grid and parameters
    load(sprintf('%s%s',Datadir,matFiles(k_t).name));
            S.w(1,1,1)=0;
        S.w(1,1,2)=0;
        S.w(1,end,1)=0;
        S.w(1,end,2)=0;
        S.w(end,1,1)=0;
        S.w(end,1,2)=0;
        S.w(end,end,1)=0;
        S.w(end,end,2)=0;
    for i =1:P.Nx
        for j = 1:P.Ny
            if((S.w(j,i,1)^2+S.w(j,i,2)^2)*(1-S.phi(j,i))^2>phi_w(1)^2+phi_w(2)^2)
               phi_w(1)= (1-S.phi(j,i))*S.w(j,i,1);
               phi_w(2)= (1-S.phi(j,i))*S.w(j,i,2);
            end
            if(S.T(j,i)>T_max)
               T_max = S.T(j,i); 
            end
            if(S.phi(j,i)>phi_max)
               phi_max = S.phi(j,i); 
            end
            if(S.phi(j,i)<phi_min)
               phi_min = S.phi(j,i); 
            end
        end
    end
end

for k_t=1:10:length(matFiles)/3*2
    disp(k_t)
        % load the grid and parameters
    load(sprintf('%s%s',Datadir,matFiles(k_t).name));

        plot(1);
        % Get x,y coordinates from h
        [x y] = visualize(h);
         x = P.l*100*x; % Nondimensional scaling in cm
         y = P.l*100*y;
        hold on;
        Nx = P.Nx;
        Ny = P.Ny;
        % Color Map
        % Temp
        %hpc = pcolor(x,y,S.T*(P.T_D-P.T_0));
        %set(hpc,'edgecolor','none','facelighting','flat','facecolor','interp');
        %caxis([0 75])
        %caxis([P.T_0-273 T_max*(P.T_D-P.T_0)+P.T_0-273])
        % Phi
        hpc = pcolor(x,y,S.phi);
        set(hpc,'edgecolor','none','facelighting','flat','facecolor','interp');
        %caxis([P.phi_0 phi_max])
        caxis([phi_min phi_max])
        
        % Porosity
        %hpc = pcolor(x,y,ones(size(S.phi))-S.phi);
        %set(hpc,'edgecolor','none','facelighting','flat','facecolor','interp');
        %caxis([1-.45 (1-P.phi_0)*1.05])
        
                tic;
        grid on;
        colorbar;
        scale = 1.25;
        axis([-scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2 -scale*P.l*100*P.Lx/2 scale*P.l*100*P.Lx/2]);
        title(sprintf('Cooking Steak T_D = %d, T_0 = %d',P.T_D, P.T_0))
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
        toc
        %Plot Velocity Vector Field
        % For plotting w
        %quiver(x,y,S.w(:,:,1)/sqrt(w(1)^2+w(2)^2),S.w(:,:,2)/sqrt(w(1)^2+w(2)^2),'AutoScale','off','Color','b','LineWidth',.5);
        % For Plotting (1-phi)*w
            S.w(1,1,1)=0;
        S.w(1,1,2)=0;
        S.w(1,end,1)=0;
        S.w(1,end,2)=0;
        S.w(end,1,1)=0;
        S.w(end,1,2)=0;
        S.w(end,end,1)=0;
        S.w(end,end,2)=0;
        quiver(x,y,(ones(size(S.phi))-S.phi).*S.w(:,:,1)/sqrt(phi_w(1)^2+phi_w(2)^2)/10,(ones(size(S.phi))-S.phi).*S.w(:,:,2)/sqrt(phi_w(1)^2+phi_w(2)^2)/10,'AutoScale','off','Color','b','LineWidth',.1);
        %quiver(x,y,(ones(size(S.phi))-S.phi).*S.w(:,:,1),(ones(size(S.phi))-S.phi).*S.w(:,:,2),'Color','b','LineWidth',.5);

        ylabel('(cm)');
        xlabel(sprintf('Normalized time = %.4f, Real Time = %.4f (sec)',t,t*P.t_0),'FontSize',12)
        
        
        hold off
        frame=getframe(fig);
        
        writeVideo(writerObj1,frame);
end

close(writerObj1);
