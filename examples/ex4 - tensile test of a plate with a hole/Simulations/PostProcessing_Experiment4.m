% Post processing specific for these simulations
function PostProcessing_Experiment4
clc
clear 
close all
%% DATA
dt = 0.5e-6;
data_dump = 4;
t_choice = [200 250 300 600 650 720]*1e-6; % 30 micro-secs
n_t = ceil(t_choice/dt/data_dump);
animationOn = false;
horizon = [5 4 3 2];
n_final = 402;
%% Cement mortar
E = 5980;%3.45e3; % [MPa]
%nu = 1/3;
nu = 0.22;
rho = 2162; % [kg/m^3]
M = E*(1-nu)./(1+nu)./(1-2*nu);
% Velocities
C1 = sqrt(M*1e6/rho);
C2 = sqrt(E*1e6./2./(1+nu)/rho);
cR = RayleighC(C1,C2,nu)
%% Load data from files
%cd SIMULATIONS
% -- LPS 
for ii = 1:length(horizon)
    filename = strcat('sim_m4_d',int2str(horizon(ii)),'LPS_holeBC_v1');
    load(filename)
    x_cell{ii} = x;
    phi_LPS{ii}= phi;
end


%% Plot crack evolution
%damagePlot(x,phi_LPS{end},n_t,dt*data_dump,animationOn)

%% Plot crack tip speed
figure
LineStyleCustom = {'-';'--';'-';'-.'};
for ii=1:length(horizon)
    [tt,V_LPS] = trackCrack(x_cell{ii},phi_LPS{ii},n_final,dt*data_dump);
    plot(tt*10^6,V_LPS/cR,'LineStyle',LineStyleCustom{ii},'LineWidth',1.5,'DisplayName',strcat('Horizon = ',num2str(horizon(ii)),'mm'))
    hold all
    set(gca,'FontSize',13);
end
%plot(tt*10^6,2/3*ones(size(tt)),'-*','LineWidth',1.5,'DisplayName','2/3 c_r')
legend
xlabel('Simulation time (\mu sec)')
ylabel('Vf/c_r')
set(gca,'FontSize',13)
grid on
hold off

%% Strain Energy Density
%energyPlot(x,energy,n_t,animationOn,dt*data_dump,'LPS')

end



function damagePlot(x,phi,n,dt,animationOn)
    if animationOn
        tt = 1:2:size(phi,2);
    else
        tt = size(phi,2);
    end
    % Scatter
    figure
    for tt = tt
    scatter(x(:,1),x(:,2),[],phi(:,tt),'filled')
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 0.4]);
    xlabel x
    ylabel y
    set(gca,'FontSize',15)
    axis equal
    grid on
    title(strcat('Damaged configuration, t = ',num2str(tt*dt)))
    pause(0.001)
    end
    %% PICTURES
    for nn = 1:length(n)
    figure
    scatter(x(:,1),x(:,2),[],phi(:,n(nn)),'filled');
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 0.4]);
    xlabel x
    ylabel y
    set(gca,'FontSize',15)
    axis equal
    %grid on
    title('Damage map')
    set(gca,'FontSize',15)
    end
end

function [tt,V_l] =  trackCrack(x,phi,n_final,dt)
    data_dump = 10; % Every 80 timesteps
    if n_final < data_dump
       data_dump = n_final; 
    end
    samples = floor(n_final/data_dump);
    V_l = zeros(samples+1,1)
    tt = (0:samples)*data_dump*dt;
    for n_samp = 0:samples
        n = 1+data_dump*n_samp;
        set_dam = find(phi(:,n) > 0.30); % Binary of greater than 0.35 damage index
        x_dam = x(set_dam,:); % Set of probable tips
        x_max = max(x_dam(:,1));
        if ~isempty(set_dam)
            damx_ind = (x_dam(:,1) == x_max);% - 1e-12);
            set_dam = set_dam(damx_ind);
            if length(set_dam) > 1
                y_max = max(x_dam(damx_ind,2));
                damy_ind = (x_dam(damx_ind,2) > y_max - 1e-12);
                set_dam = set_dam(damy_ind);
                tip(2) = set_dam;
            else
                tip(2) = set_dam;
            end
            if n_samp > 1 && (norm(x(tip(2),:)-x(tip(1),:)) < 0.010 || V_l(n_samp) > 0)
                V_l(n_samp+1) = norm(x(tip(2),1)-x(tip(1),1))/dt/data_dump;
                tip(1) = tip(2);
            else
                tip(1) = set_dam;
            end
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
        animation = 1:2:size(energy.W,2);
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
    figure
    hh = surf(X,Y,WW(:,:,n(nn)));
    xlabel x
    ylabel y
    axis equal
    title('Strain Energy Density')
    set(gca,'FontSize',15)
    set(hh,'EdgeColor', 'none');
    view(2)
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

function CR = RayleighC(C1,C2,nu)
    CR = C1;
    for ii = 1:length(C1)
        fun = @(CR) (2-CR^2/C2(ii)^2)^2 - 4*sqrt((1-CR^2/C1(ii)^2)*(1-CR^2/C2(ii)^2));
        x0 = C2(ii);
        CR_comp = fsolve(fun,x0);
        CR(ii) = real(CR_comp);
    end
    %cc = sqrt((30.876 - 14.876*nu - sqrt(224.545376*nu^2 - 93.122752*nu + 124.577376))/26/(1-nu));
    %CR = cc*C2;
end
