% Post processing specific for these simulations
function PostProcessing_Experiment2
clc
clear 
close all
%% DATA
dt = 0.02e-6;
n_final = 1000; % 30 micro-secs
animationOn = false;
%% Load data from files
%cd SIMULATIONS
% -- PMB DTT
cd 'PMB DTT'
load('sim_m4_d1_j1') % 2 MPA
%load('sim_m4_d1_j1') % 4 MPA
phi_PMB_DTT = phi;
% -- PMB 
cd ../PMB
%load('sim_m4_d1_j1_tt2') % 2 MPA
load('sim_m4_d1_j1_tt4') % 4 MPA
phi_PMB = phi;
% -- LPS 
cd ../LPS
load('sim_m4_d1LPS_2')
[phi_LPS]= phi;
% -- Lipton
cd ../Lipton
load('sim_m4_d1LSJT')
[phi_Lipton]= phi;


%% Plot crack evolution
damagePlot(x,phi_PMB,'PMB',dt,animationOn,n_final)
damagePlot(x,phi_PMB_DTT,'PMB DTT',dt,animationOn,n_final)
damagePlot(x,phi_LPS,'LPS',dt,animationOn,n_final)
damagePlot(x,phi_Lipton,'Lipton',dt,animationOn,n_final)
%% Plot crack velocity
[tt0,V_PMB] = trackCrack(x,phi_PMB,n_final,dt);
[tt1,V_PMB_DTT] = trackCrack(x,phi_PMB_DTT,n_final,dt);
[tt2,V_LPS] = trackCrack(x,phi_LPS,n_final,dt);
[tt3,V_Lipton] = trackCrack(x,phi_Lipton,n_final,dt);
cR_1 = 3102; % m/s
cR_2 = 3180; % m/s
% Plot the velocity
figure
plot(tt0*10^6,V_PMB/cR_2,'--','LineWidth',1.5)
hold on
plot(tt1*10^6,V_PMB_DTT/cR_2,'-d','LineWidth',1.5)
plot(tt2*10^6,V_LPS/cR_2,'-s','LineWidth',1.5)
plot(tt3*10^6,V_Lipton/cR_2,'-*','LineWidth',1.5)
plot(tt1*10^6,ones(size(tt1))*2/3,'-','LineWidth',1.5)
legend('PMB','PMB DTT','LPS','Lipton','2/3 C_r')
%legend('PMB DTT','LPS','Lipton','Experimental tests')
xlabel('Simulation time (\mu sec)')
ylabel('Vf/c_r')
set(gca,'FontSize',13)
grid on
hold off

cd ../

l = 0.01;
theta_LPS = angleEvaluation(x,phi_LPS,l)
theta_Lipton = angleEvaluation(x,phi_Lipton,l)
theta_PMB = angleEvaluation(x,phi_PMB,l)
theta_PMBDTT = angleEvaluation(x,phi_PMB_DTT,l)


end

function [tt,V_l] =  trackCrack(x,phi,n_final,dt)
    data_dump = 50; % Every 80 timesteps
    if n_final < data_dump
       data_dump = n_final; 
    end
    samples = floor(n_final/data_dump);
    V_l = zeros(samples+1,1);
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
            if n_samp > 1
                V_l(n_samp+1) = norm(x(tip(2),1)-x(tip(1),1))/dt/data_dump;
                tip(1) = tip(2);
            else
                tip(1) = set_dam;
            end
        end           
    end
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
        animation = 1:5:size(phi,2);
    else
        animation = n_final;
    end
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
    pause(0.001)
    end
end

function theta = angleEvaluation(x,phi,l)
h = norm(x(1,:) - x(2,:));
n_final = size(phi,2);
x_initial = max(x(phi(abs(x(:,2))<=h/2,n_final)>0.45,1));
y_initial = 0;
x_tip = [x_initial x_initial+l];
y_tip = [y_initial];
column_index = find(abs(x(:,1)-x_tip(2)) < h/2 & x(:,2)>0); %checkFamily(x,[],column_index,true)
index_tip  = column_index(phi(column_index,n_final) == max(phi(column_index,n_final)));
y_tip = [y_tip, x(index_tip,2)];

% Figure
damagePlot(x,phi,"angle",0,false,size(phi,2))
hold on
plot(x_tip,y_tip,'m','LineWidth',2)

% Find angle
theta = atan(diff(y_tip)/diff(x_tip))*180/pi;


end