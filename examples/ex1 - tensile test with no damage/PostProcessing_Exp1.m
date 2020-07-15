% Post processing
function PostProcessing_Exp1()
close all
%% DATA
horizon = [1]*1e-3; % [m]
E = 72e3; % [MPa]
nu = 1/3;%0.22;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
M = E*(1-nu)/(1+nu)/(1-2*nu);
dt = 0.02e-6; % 0.02 micro-secs
% Velocities
C1_1 = sqrt(M*1e6/rho);
C2_1 = sqrt(E*1e6/2/(1+nu)/rho);
nu = 0.22;
M = E*(1-nu)/(1+nu)/(1-2*nu);
C1_2 = sqrt(M*1e6/rho);
C2_2 = sqrt(E*1e6/2/(1+nu)/rho);
CR_1 = RayleighC(C1_1,C2_1,nu);
CR_2 = RayleighC(C1_2,C2_2,nu);
%% LOAD STAGE
cd 'Simulations/PMB'
horizon = [1 2 3 4];
%x_n = cell(size(horizon));
for ii = 1:length(horizon)
    filename = strcat('sim_m4_d',num2str(horizon(ii)),'PMB DTT.mat');
    load(filename)
    energy_PMB(ii) = {energy}; x_PMB(ii) = {x};
end

cd '../LPS 2D'
for ii = 1:length(horizon)
    filename = strcat('sim_m4_d',num2str(horizon(ii)),'LPS 2D.mat');
    load(filename)
    energy_LPS(ii) = {energy}; x_LPS(ii) = {x};
end

cd '../Lipton'
for ii = 1:length(horizon)
    filename = strcat('sim_m4_d',num2str(horizon(ii)),'Lipton Free Damage.mat');
    load(filename)
    energy_Lipton(ii) = {energy}; x_Lipton(ii) = {x};
end

cd ../../
energyPlot(x_PMB{1},energy_PMB{1},1501);
energyPlot(x_LPS{1},energy_LPS{1},1501);
energyPlot(x_Lipton{1},energy_Lipton{1},1501);

%% 
%u = threeDModification(x,u_n,idb);
for ii = 1:length(horizon)
    n_PMB(ii) = energyPlot(x_PMB{ii},energy_PMB{ii},1501);
    n_LPS(ii) = energyPlot(x_LPS{ii},energy_LPS{ii},1501);
    n_Lipton(ii) = energyPlot(x_Lipton{ii},energy_Lipton{ii},1501);
end
close all
n_PMB = [358 n_PMB]; n_LPS = [363 n_LPS]; n_Lipton = [367 n_Lipton];
n_Abaqus = [361 353];
% Manually adding the 0.5 result
t_PMB = dt; t_LPS = n_LPS*dt; t_Lipton = n_Lipton*dt;
horizon_PMB = [0.5 horizon];

dL = 0.02; % [m]

V_PMB = dL./t_PMB; V_LPS = dL./t_LPS; V_Lipton = dL./t_Lipton;

% figure 
% plot(horizon_PMB,V_PMB,'-.','DisplayName','PMB')
% hold on
% plot(horizon,V_LPS,'-*','DisplayName','LPS')
% plot(horizon,V_Lipton,'-s','DisplayName','LJS')
% plot(horizon,CR_1*ones(size(horizon)),'--','DisplayName','CR \nu = 1/3')
% plot(horizon,CR_2*ones(size(horizon)),'-','DisplayName','CR \nu = 0.22')
% xlabel('Horizon (mm)')
% ylabel('Rayleigh speeds (m/sec)')
% grid on
% legend

% Figure - Frequency
w_pmb = 1./n_PMB/dt/2; 
w_lps = 1./n_LPS/dt/2;
w_lipton = 1./n_Lipton/dt/2;
w_Abaqus = 1./n_Abaqus/dt/2;
horizon = [0.5 horizon];
figure 
plot(horizon,w_pmb,'-.','DisplayName','PMB','LineWidth',1.5)
hold on
plot(horizon,w_lps,'-*','DisplayName','LPS','LineWidth',1.5)
plot(horizon,w_lipton,'-s','DisplayName','LJS','LineWidth',1.5)
plot(horizon,w_Abaqus(1)*ones(size(horizon)),'--','DisplayName','Abaqus \nu = 0.22','LineWidth',1.5)
plot(horizon,w_Abaqus(2)*ones(size(horizon)),'-','DisplayName','Abaqus \nu = 0.3','LineWidth',1.5)
xlabel('Horizon (mm)')
ylabel('Strain energy frequency (Hz)')
grid on
legend

error = [abs(w_pmb(1)-w_Abaqus(2))/w_Abaqus(2) abs(w_lps(1)-w_Abaqus(1))/w_Abaqus(1) abs(w_lipton(1)-w_Abaqus(1))/w_Abaqus(1)];
end

function CR = RayleighC(C1,C2,nu)
    fun = @(CR) (2-CR^2/C2^2)^2 - 4*sqrt((1-CR^2/C1^2)*(1-CR^2/C2^2));
    x0 = C2;
    CR = fsolve(fun,x0);
    CR = real(CR);
    %cc = sqrt((30.876 - 14.876*nu - sqrt(224.545376*nu^2 - 93.122752*nu + 124.577376))/26/(1-nu));
    %CR = cc*C2;
end

function loc = energyPlot(x,energy,n_final)
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
    title('Strain Energy Density')
    set(gca,'FontSize',15)
   %% Total energy evolution
   int = energy.W + energy.KE - energy.EW; % Total energy for each point
   t = 1:1:n_final; t = t*0.02e-16;% Number of time steps
   [pks,loc] = findpeaks(-sum(energy.KE(:,1:n_final)));
   figure
   plot(t,sum(energy.W(:,1:n_final)),'--','LineWidth',1.5)
   hold on
   plot(t,sum(energy.KE(:,1:n_final)),'-.','LineWidth',1.5)
   plot(t,sum(energy.EW(:,1:n_final)),'-','LineWidth',1.5)
   plot(t,sum(int(:,1:n_final)),'-','LineWidth',2.0)
   %plot(loc,-pks,'or')
   legend('Strain energy','Kinetic energy','External work','Total internal energy','Location','bestoutside')
   xlabel('Time step n')
   ylabel('Energy (J/m)')
   set(gca,'FontSize',15)
   grid on

   %
   loc = loc(1);
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