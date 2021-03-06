Datadir ='.\perm1e-17_diff_CNM_thick\'
infilebase =  sprintf('%s/steak_default',Datadir);

    % load grid and parameters
    %
    paramfile = sprintf('%s.param.mat',infilebase);
    load(paramfile);
    
    filePattern = fullfile(Datadir, 'steak_default.*.mat');
matFiles = dir(filePattern);
names = matFiles.name;
numFilesTenMin = floor(1*60/P.t_0/P.dt/P.outevery);
hold on;
counter = 1;
for k_t=1:numFilesTenMin:length(matFiles)
    disp(k_t)
        % load the grid and parameters
    load(sprintf('%s%s',Datadir,matFiles(k_t).name));
            [x y] = visualize(h);
         x = P.l*100*x; % Nondimensional scaling in cm
         y = P.l*100*y;
         XX(counter,:,:)=x;
         YY(counter,:,:)=y;
             T_all(counter)=S.T(end/2,end/2+.5);
             Time(counter)=t;
    % Horizontal Slice, Moisture Conent  
    %plot(x(end/2,:),(ones(size(S.phi(end/2,:)))-S.phi(end/2,:))./(ones(size(S.phi(end/2,:)))+.3*S.phi(end/2,:)),'LineWidth',2);
    % Vertical Slice, Moisture Conent  
    %plot(y(:,ceil(end/2)),100*(ones(size(S.phi(:,ceil(end/2))))-S.phi(:,ceil(end/2)))./(ones(size(S.phi(:,ceil(end/2))))+.3*S.phi(:,ceil(end/2))),'LineWidth',2);
    %    if(counter == 4)
    %plot(y(:,ceil(end/2)),(ones(size(S.phi(:,ceil(end/2))))-S.phi(:,ceil(end/2)))./(ones(size(S.phi(:,ceil(end/2))))+.3*S.phi(:,ceil(end/2))),'--','Color',	[0, 0.4470, 0.7410],'LineWidth',2);
    %    end
    %    if(counter ==7)
    %plot(y(:,ceil(end/2)),(ones(size(S.phi(:,ceil(end/2))))-S.phi(:,ceil(end/2)))./(ones(size(S.phi(:,ceil(end/2))))+.3*S.phi(:,ceil(end/2))),'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
    %    end
    counter = counter+1
    %Moisture Content Orig B.C. + Diffusion, Horizontal, T_D = 175 \circ C, c_0 Bengston
    %Moisture Content Orig B.C. + Diffusion, Vertical, T_D = 175 \circ C, c_0 Bengston
    % Temperature 
    %plot(x(end/2,:),S.T(end/2,:)*(P.T_D-P.T_0)+(P.T_0-273)*ones(size(S.T(end/2,:))),'LineWidth',2)
    %plot(y(:,ceil(end/2)),S.T(:,ceil(end/2))*(P.T_D-P.T_0)+(P.T_0-273)*ones(size(S.T(:,ceil(end/2)))),'LineWidth',2)
    %if(counter == 3)
    %plot(y(:,ceil(end/2)),S.T(:,ceil(end/2))*(P.T_D-P.T_0)+(P.T_0-273)*ones(size(S.T(:,ceil(end/2)))),'--','Color',	[0, 0.4470, 0.7410],'LineWidth',2)
    %end
    %if(counter == 4)
    %plot(y(:,ceil(end/2)),S.T(:,ceil(end/2))*(P.T_D-P.T_0)+(P.T_0-273)*ones(size(S.T(:,ceil(end/2)))),'--','Color',	[0.8500, 0.3250, 0.0980],'LineWidth',2)
    %end
    %if(counter == 5)
    %plot(y(:,ceil(end/2)),S.T(:,ceil(end/2))*(P.T_D-P.T_0)+(P.T_0-273)*ones(size(S.T(:,ceil(end/2)))),'--','Color',	[0.9290, 0.6940, 0.1250],'LineWidth',2)
    %end
    % pi_sw
    %    PI_sw = P.mu_0/P.t_0*(pi_el(P,S.T,S.phi) + pi_mix(P,S.T,S.phi));
    %plot(x(end/2,:),PI_sw(end/2,:),'LineWidth',2)
    %plot(y(:,17),PI_sw(:,17),'LineWidth',2)

end
return;
% Moisture Content DATA
if(0)
   x = [4 8 12 16 20 26 32 40 47.5 50 55];
    x=x/10-2.75*ones(size(x));
    x=x*YY(3,end,ceil(end/2))/2.75;
    n = [48 51 56 66 72 71 70 66 62 56 52];
    n=n/100;
    plot(x,n,'o','Color',[0, 0, 1],'LineWidth',3); 
   x = [4 8 16 20 32 38 42 47.5 50 55];
    x=x/10-2.75*ones(size(x));
    x=x*YY(3,end,ceil(end/2))/2.75;
    n = [40 45 47 53 52 51 50.5 48 43 41];
    n=n/100;
    plot(x,n,'o','Color',[1, 0, 0],'LineWidth',3); 
end

if(0)
x = [4 8 12 16 20 26 32 38 40 42.5 45 47.5 50 52.5 55];
x=x/10-2.75*ones(size(x));
x=x*YY(1,end,ceil(end/2))/2.75;
T = [42 36 21 16 11 9.5 11 16 20 24 31 40 44 58 60];
plot(x,T,'o','Color',[0, 0, 1],'LineWidth',3);
x = [2 5 8 12 15 18 20 23 26.5 28 29.5 32.5 35 40 44 48 52 54 55];
x=x/10-2.75*ones(size(x));
x=x*YY(2,end,ceil(end/2))/2.75;
T = [80 61 51 43 39 29 24 23 21 21 23 24 29 41 50 58 62 67 80];
plot(x,T,'o','Color',[1, 0, 0],'LineWidth',3);
x = [3 6 10 14 18 26 30 36 39 42 46 52 55];
x=x/10-2.75*ones(size(x));
x=x*YY(3,end,ceil(end/2))/2.75;
T = [78 70 61 56 46 41 39 41.5 50 58 66 75 83];
plot(x,T,'o','Color',[0.75, 0.75, 0],'LineWidth',3);
end
%Temperature T_D = 175 \circ C \kappa = 1e-16 Diff = 7e-9
grid on;
xlabel('(cm)')
%ylabel('\circ C')
title('Moisture content, vertical slice')
ylabel('% mass')
%legend('0','10','20','30','40','50','60','70','80','90','100','110','120');