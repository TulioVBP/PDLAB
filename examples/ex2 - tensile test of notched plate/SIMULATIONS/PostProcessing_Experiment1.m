% Post processing specific for these simulations
function PostProcessing_Experiment1
clc
clear 
close all
%% DATA
dt = 0.02e-6;
data_dump = 5;
t_choice = [0.001 8 10 12 14 16 18 21]*1e-6; % 30 micro-secs
n_t = ceil(t_choice/dt/data_dump);
animationOn = false;
%% Load data from files
%cd SIMULATIONS
% -- LPS 
cd Lipton
load('sim_m4_d1Lipton_data_dump')
[phi_LPS]= phi;


%% Plot crack evolution
damagePlot(x,phi_LPS,'Lipton',dt*data_dump,animationOn,n_t)

%% Strain Energy Density
%energyPlot(x,energy,n_t,animationOn,dt*data_dump,'LPS')

cd ../
end



function damagePlot(x,phi,modelname,dt,animationOn,n_final)
    b = max(x(:,1)) - min(x(:,1));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
    %for n = 1:size(phi,2)
    PHI = zeros([size(X),size(phi,2)]);
    sz = flip(size(X));
    for n = 1:size(phi,2)
         PHI(:,:,n)= reshape(phi(:,n),sz)';
    end
    if animationOn
        animation = 1:2:size(phi,2);
        figure
        for n = animation
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
        pause(0.01)
        end
    else
        animation = n_final;
        for n = animation
        figure
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
        pause(0.01)
        end
    end
end

function energyPlot(x,energy,n,animationOn,dt,modelname)
   %% Strain energy density
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
    WW = zeros([size(X),size(energy.W,2)]);
    sz = flip(size(X));
    for t = 1:size(energy.W,2)
        WW(:,:,t)= reshape(energy.W(:,t),sz)';
    end
    %% ANIMATION
    if animationOn
        animation = 1:5:size(energy.W,2);
    else
        animation = n_final;
    end
    figure
    for nA = animation
    hh = pcolor(X,Y,WW(:,:,nA));
    axis equal
    set(hh,'EdgeColor', 'none');
    xlabel x
    ylabel y
    title(strcat('Strain energy density - ',modelname,'- t = ',num2str(nA*dt),' secs'))
    set(gca,'FontSize',15)
    c = jet(1000);
    colormap(c);
    colorbar
    %caxis([0 0.5]);
    pause(0.01)
    end
    %% PICTURES
    for nn = 1:length(n)
    avg_value = sum(energy.W(:,n(nn)))/length(energy.W(:,n(nn)));
    n_down = sum(energy.W(:,n(nn))<1.25*avg_value);
    n_up = length(energy.W(:,n(nn))) - n_down;
    top_color = (n_down*avg_value + n_up*max(energy.W(:,n(nn))))/(n_down+n_up);
    figure
    hh = surf(X,Y,WW(:,:,n(nn)));
    xlabel x
    ylabel y
    axis equal
    title('Strain Energy Density')
    set(gca,'FontSize',15)
    set(hh,'EdgeColor', 'none');
    view(2)
    colorbar
    caxis([0 top_color]);
    end

   %% Total energy evolution
   int = energy.W + energy.KE - energy.EW; % Total energy for each point
   t = 1:1:n(nn); % Number of time steps
   figure
   plot(t,sum(energy.W(:,1:n(end))),'--','LineWidth',1.5)
   hold on
   plot(t,sum(energy.KE(:,1:n(end))),'-.','LineWidth',1.5)
   plot(t,sum(energy.EW(:,1:n(end))),'-','LineWidth',1.5)
   plot(t,sum(int(:,1:n(nn))),'-','LineWidth',2.0)
   legend('Strain energy','Kinetic energy','External work','Total internal energy')
   xlabel('Time step n')
   ylabel('Energy (J/m)')
   set(gca,'FontSize',15)
   grid on
   
end