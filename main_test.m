% MAIN SCRIPT TO RUN A 2D PERIDYNAMIC BOND BASED MODEL - TEST
clear all
clc
close all

%% PARAMETERS
global a b alpha c1 horizon omega Sc S0 S1 damageOn rho
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
sigma = 4; % [MPa]
dt = 0.02e-6; % [sec]
m_vec = 2; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; % Influence function option
notch_length = 0.05; % 5 cm
crackIn = [0 0.02; notch_length 0.02]; % Coordinates of the crack initial segment
damageOn = false; % True if applying damage to the model, false if not
model = "Linearized LPS bond-based";
solver = "Quasi-Static"; % "Quasi-Static"
switch model
    case "PMB"
        alpha = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        %c1 = 24*(E*1e6)/pi/horizon^3/(1-nu); % Pa/m^3 = N / m^5
        m = weightedVolume(horizon,omega);
        c1 = 6*E*1e6/m;
        T = @interactionForce_PMB;
    case "Linearized LPS bond-based"
        alpha = 0; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        m = weightedVolume(horizon,omega); 
        c1 = 6*E*1e6/m;
        T = @interactionForce_LLPSBB;
    otherwise
        disp("Chosen model is not implemented or it was mistyped");
        quit
end

%% SIMULATION
for s_index = 1:length(sigma)
    %% STRESS LOOP
    stresses = [0 sigma(s_index) 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for m_index = 1:length(m_vec)
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
        A = h^2; % Elements' area
        [ndof,idb,bc_set,bodyForce] = boundaryCondition(x,stresses,m,h,A);
        % ###### GENERATE FAMILY #######
        [family,partialAreas,maxNeigh] = generateFamily(x,horizon,m,m_index,true); % True for test
        %% SOLVER
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t = 0:dt:40e-6; % 40 micro-secs simulation
                [u_n,phi] = solver_DynamicExplicit(x,maxNeigh,t,idb,bodyForce,family,partialAreas,T,crackIn);
            case "Quasi-Static"
                n_tot = 4;
                [un,r] = solver_QuasiStatic(x,n_tot,bodyForce,idb,family,partialAreas,T,ndof,A);
            otherwise
                disp("ERROR: Solver not yet implemented.")
                quit
        end    
    end
end
PostProcessing(x,u_n,n,phi,idb);
