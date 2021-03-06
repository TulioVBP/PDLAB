% MAIN SCRIPT TO RUN A 2D PERIDYNAMIC BOND BASED MODEL
clear all
clc
close all

%% PARAMETERS
global alpha c1 horizon omega Sc S0 S1
% Material
horizon = 1e-3; % [m]
E = 72e3; % [MPa]
nu = 1/3;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
S0 = [-0.98 0.95*Sc]; % S0- and S0+
S1 = [-0.99 1.05*Sc]; % S1- and S1+
% Simulation
sigma = [0.2 2 4]; % Final simulation [MPa]
dt = 0.02e-6; % [sec]
m_vec = [2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
alpha = 0;
omega = 3; % Influence function option
c1 = 24*E/pi/horizon^3/(1-nu);


%% SIMULATION
for s_index = 1:1%length(sigma)
    %% STRESS LOOP
    stresses = [0 sigma(s_index) 0]; % [sigma_x, sigma_y, tau_xy]
    for m_index = 1:1%length(m+vec)
        %% HORIZON NUMBER LOOP
        % ###### GENERATE MESH #######
        h = h_vec(m_index); % grid spacing [m]
        m = m_vec(m_index); 
        a = 0.04; % height [m]
        b = 0.10; % length [m]
        N = floor(a/h+1/2) + 1; % Number of rows
        M = floor(b/h+1/2) + 1; % Number of collumns
        x = zeros(N*M,2);
        for ii = 1:N
            for jj = 1:M
                x((ii-1)*M+jj,:) = [(jj-1)*h,(ii-1)*h];
            end
        end
        A = h^2*ones(length(x),1); % Elements' area
        % Plot the mesh
%         figure(1)
%         plot(x(:,1),x(:,2),'o')
%         grid on
%         xlabel('x[m]')
%         ylabel('y[m]')
        % ###### GENERATE FAMILY #######
        [family,partialAreas,maxNeigh] = generateFamily(x,horizon,m,m_index);
        %% REAL TIME SIMULATION
        dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
        t = 0:dt:150e-6; % 1 secs simulation
        % ########### INITIALIZE SIMULATION MATRICES ##############
        M = rho * eye(length(x)); % Inertial matrix
        Minv = 1/rho * eye(length(x)); % Inverse of the inertial matrix
        S_max = zeros(length(x),maxNeigh);
        phi = ones(length(x),length(t)); % No initial damage
        % Initial condition
        u_n = zeros(length(x),2*length(t));
        %u = zeros(size(x));
        v_n = zeros(length(x),2*3); % Velocity Verlet matrix [i i+1/2 i+1]: by initializing it to zero, the rigid motion is eliminated.ffz
        bn = bodyForce(x,stresses,m,h,A);
        for n = 1:length(t)-1
            %% TIME LOOP
            fn = zeros(size(x));
            for ii = 1:length(x)
                % Loop on the nodes
               for neig_index=1:size(family,2)
                   % Loop on their neighbourhood
                   jj = family(ii,neig_index); 
                   if jj == 0 % Check if the jj neighbor is a valid node
                       break
                   elseif ii~=jj % Check if jj is a different node (shouldn't be necessary as family already accounts for it)
                       [fij,S_max(ii,neig_index)] = interactionForce(x(ii,:),x(jj,:),u_n(ii,(2*n-1):(2*n)),u_n(jj,(2*n-1):(2*n)),S_max(ii,neig_index));
                       Vj = partialAreas(ii,neig_index);
                       fn(ii,:) = fn(ii,:) + fij*Vj;
                   end
               end
               phi(ii,n) = damageIndex(x,u_n(:,(2*n-1):(2*n)),family(ii,:),partialAreas(ii,:),ii); % Damage index
            end
            % ############ VELOCITY VERLET ALGORITHM ###############
            %a_n = dt/2*M\(fn + bn); % Acceleration
            % Step 1 - Midway velocity
            v_n(:,3:4) = v_n(:,1:2) + dt/2*Minv*(fn + bn); % V(n+1/2)
            u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4);
            v_n(:,5:6) = v_n(:,3:4) + dt/2*Minv*(fn + bn); % V(n+1) - unnecessary
            v_n(:,1:2) = v_n(:,5:6); % For the next iterative loop
            % ############ COUNTING THE PROCESSING TIME #############
            disp("Stress: "+ int2str(s_index) + " - Mesh: " +int2str(m_index) + " - Percentage of the process: " + num2str(n/length(t)*100) + "%")
            %pause
        end
    end
end
