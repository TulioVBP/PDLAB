% Post processing specific for these simulations
function PostProcessing_Experiment2
clc
clear 
close all
%% DATA
dt = 0.5e-6;
n_final = 2001; % 20 micro-secs
%% Load data from files
%cd SIMULATIONS
% PMB DTT
cd 'PMB DTT'
load('sim_m4_d20_j1_traction')
phi_PMB_DTT = phi;
x_PMB = x;
u_PMB = u_n;
idb_PMB = idb;
energy_PMB = energy;
cd ..

cd 'LPS'
load('sim_m4_d20_j1_traction')
phi_LPS = phi;
u_LPS = u_n;
x_LPS = x;
idb_LPS = idb;
energy_LPS = energy;
cd ..

%% Plot crack evolution
%damagePlot(x,phi_PMB_DTT,'PMB DTT',dt)

%% Plot deformed shape
u_PMB = threeDModification(x_PMB,u_PMB,idb_PMB);
u_LPS = threeDModification(x_LPS,u_LPS,idb_LPS);

deformedShape(x_PMB,u_PMB,phi_PMB_DTT,'PMB DTT',dt,true)
deformedShape(x_LPS,u_LPS,phi_LPS,'LPS',dt,true)

%% Energy plot
energyPlot(x_PMB,energy_PMB,n_final)
energyPlot(x_LPS,energy_LPS,n_final)

end

function damagePlot(x,phi,modelname,dt)
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-12, min(x(:,2)):h:max(x(:,2))+1e-12);
    %for n = 1:size(phi,2)
    PHI = zeros([size(X),size(phi,2)]);
    sz = flip(size(X));
    for n = 1:size(phi,2)
         PHI(:,:,n)= reshape(phi(:,n),sz)';
    end
    
    figure
    for n = 1:5:size(phi,2)
    hh = pcolor(X,Y,PHI(:,:,n));
    axis equal
    set(hh,'EdgeColor', 'none');
    xlabel x
    ylabel y
    title(strcat('Damage index -',modelname,'- t = ',num2str(n*dt),' secs'))
    set(gca,'FontSize',15)
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 0.5]);
    pause(0.001)
    end
end

function deformedShape(x,u,phi,modelname,dt,animationOn)
    b = max(x(:,1)) - min(x(:,1));
    a = max(x(:,2)) - min(x(:,2));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
    V = zeros([size(X) size(u,2)]);
    W = zeros([size(X) size(u,2)]);
    sz = flip(size(X));
    for n = 1:size(u,3)
        V(:,:,n) = reshape(u(:,1,n),sz)';
        W(:,:,n) = reshape(u(:,2,n),sz)';
    end
    %% Plot displacement in the mesh
    scaleFactor = 1e3;
    % Mesh
    figure
    if animationOn
        for n = 1:5:size(u,3)
            %mesh(X+scaleFactor*V(:,:,n),Y+scaleFactor*W(:,:,n),'LineWidth',1)
            scatter(x(:,1)+scaleFactor*u(:,1,n),x(:,2)+scaleFactor*u(:,2,n),[],phi(:,n),'filled')
            hidden off
            xlabel x
            ylabel y
            set(gca,'FontSize',15)
            title(strcat('Damage index -',modelname,'- t = ',num2str(n*dt),' secs'))
            c = jet(1000);
            colormap(c);
            colorbar
            caxis([0 0.4]);
            axis equal
            view(2)
            pause(0.1)
        end
    else
        scatter(x(:,1)+scaleFactor*u(:,1,size(u,3)),x(:,2)+scaleFactor*u(:,2,size(u,3)),[],phi(:,size(u,3)),'filled')
        hidden off
        xlabel x
        ylabel y
        set(gca,'FontSize',15)
        title(strcat('Damage index -',modelname,'- t = ',num2str(n*dt),' secs'))
        c = jet(1000);
        colormap(c);
        colorbar
        caxis([0 0.4]);
        axis equal
        view(2)
    end
end

function u = threeDModification(x,u_2D,idb)
    if size(u_2D,1) == 2*length(x)
        u = zeros(length(x),2,size(u_2D,2));
        ii = (1:length(x))';
            dofi = [idb(2*ii-1) idb(2*ii)];
            for n = 1:size(u_2D,2)
                u_temp = u_2D(:,n);
                u(:,:,n) = u_temp(dofi);
            end
    else
        u = u_2D;
    end
    disp('Displacement vector converted');
end

function energyPlot(x,energy,n_final)
   %% Strain energy density
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
    WW = zeros(size(X));
    for jjj = 1:size(X,2)
        for iii = 1:size(Y,1)
            %ind = find(x(:,1) == X(1,jj) & x(:,2) == Y(ii,1)); not
            %reliable
            coord = [X(1,jjj) Y(iii,1)];
            for ii = 1:length(x)
               if x(ii,1) < coord(1) + 1e-10 && x(ii,1) > coord(1) - 1e-6  && x(ii,2) < coord(2) + 1e-10 && x(ii,2) > coord(2) - 1e-10
                   ind = ii;
                   break;
               end
            end
            if ~isempty(ind)
                WW(iii,jjj) = energy.W(ind,n_final);
            else
                disp('Error: problem with the matching condition')
            end
        end        
    end
    figure
    surf(X,Y,WW)
    xlabel x
    ylabel y
    axis equal
    title('Strain Energy Density')
    set(gca,'FontSize',15)
   %% Total energy evolution
   int = energy.W + energy.KE - energy.EW; % Total energy for each point
   t = 1:1:n_final; % Number of time steps
   figure
   plot(t,sum(energy.W(:,1:n_final)),'--','LineWidth',1.5)
   hold on
   plot(t,sum(energy.KE(:,1:n_final)),'-.','LineWidth',1.5)
   plot(t,sum(energy.EW(:,1:n_final)),'-','LineWidth',1.5)
   plot(t,sum(int(:,1:n_final)),'-','LineWidth',2.0)
   legend('Strain energy','Kinetic energy','External work','Total internal energy')
   xlabel('Time step n')
   ylabel('Energy (J/m)')
   set(gca,'FontSize',15)
   grid on

   
end