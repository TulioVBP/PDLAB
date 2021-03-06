%% Script to run alternative test on the implementation
clear all
clc
close all

%% PARAMETERS
global a b alpha c1 horizon omega Sc S0 S1 damageOn rho crackIn model
% Material
horizon = 5e-3; % [m]
E = 72e3; % [MPa]
nu = 1/3;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
S0 = [-0.98 0.95*Sc]; % S0- and S0+
S1 = [-0.99 1.05*Sc]; % S1- and S1+
% Simulation
sigma = 4; % [MPa]
dt = 0.02e-6; % [sec]
m_vec = [2 4]; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; % Influence function option
notch_length = 0.05; % 5 cm
crackIn = [0 0.02; notch_length 0.02]; % Coordinates of the crack initial segment
damageOn = false; % True if applying damage to the model, false if not
model.name = "Linearized LPS bond-based"; % "PMB", "Linearized LPS bond-based"
solver = "Quasi-Static"; % "Quasi-Static", "Dynamic/Explicit"
switch model.name
    case "PMB"
        alpha = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        %c1 = 24*(E*1e6)/pi/horizon^3/(1-nu); % Pa/m^3 = N / m^5
        m = weightedVolume(horizon,omega);
        c1 = 6*E*1e6/m;
        T = @interactionForce_PMB;
        model.linearity = false;
    case "Linearized LPS bond-based"
        alpha = 0; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        m = weightedVolume(horizon,omega); 
        c1 = 6*E*1e6/m;
        T = @interactionForce_LLPSBB;
        model.linearity = false; %true;
    otherwise
        disp("Chosen model is not implemented or it was mistyped");
        pause
end

%% SIMULATION
for s_index = 1:length(sigma)
    % STRESS LOOP
    stresses = [0 sigma(s_index) 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for m_index = 1%:length(m_vec)
        % HORIZON NUMBER LOOP
        %% ###### GENERATE MESH #######
        h = h_vec(m_index); % grid spacing [m]
        m = m_vec(m_index); 
        a = 0.04/2; % height [m]
        b = 0.10/2; % length [m]
        N = floor(a/h+1/2) + 1; % Number of rows
        M = floor(b/h+1/2) + 1; % Number of collumns
        x = zeros(N*M,2);
        for ii = 1:N
            for jj = 1:M
                x((ii-1)*M+jj,:) = [(jj-1)*h,(ii-1)*h];
            end
        end
        A = h^2; % Elements' area
        %% Boundary conditions
        [ndof,idb,bc_set,bodyForce,noFailZone] = boundaryCondition(x,stresses,m,h,A);
        %% ###### GENERATE FAMILY #######
        [family,partialAreas,maxNeigh] = generateFamily(x,horizon,m,m_index,true); % True for test
        % SOLVER
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t = 0:dt:40e-6; % 40 micro-secs simulation
                n_tot = length(t);
                [u_n,phi] = solver_DynamicExplicit(x,maxNeigh,t,idb,bodyForce,bc_set,family,partialAreas,T,crackIn,noFailZone);
            case "Quasi-Static"
                for opt = 1:2
                    n_tot = 2;
                    [u_n(:,:,opt),r(opt)] = solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,T,ndof,A,opt);
                end
            otherwise
                disp("ERROR: Solver not yet implemented.")
                pause
        end    
    end
end

%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        PostProcessing(x,u_n,n_tot,idb,phi);
    case "Quasi-Static"
        for opt = 1:2
            PostProcessing(x,u_n(:,:,opt),n_tot,idb);
        end
        % Evaluate the error
        error = norm(u_n(:,n_tot,1) - u_n(:,n_tot,2));
    otherwise
end
