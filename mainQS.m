% MAIN SCRIPT TO RUN A 2D PERIDYNAMIC BOND BASED MODEL - TEST
clear all
clc
close all

%% PARAMETERS
global c1 c2 nu Sc S0 S1 damageOn rho crackIn model
% --- Material
horizon = 1e-3; % [m]
E = 72e3; % [MPa]
nu = 1/3;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
S0 = [-0.98 0.95*Sc]; % S0- and S0+
S1 = [-0.99 1.05*Sc]; % S1- and S1+
% ---- Simulation
sigma = 4; % [MPa]
m_vec = [2 4]; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; alfa = 1;% Influence function options
par_omega = [horizon omega alfa];
notch_length = 0.05; % 5 cm
crackIn = [0 0.02; notch_length 0.02]; % Coordinates of the crack initial segment
% ---- MODEL
damageOn = false; % True if applying damage to the model, false if not
model.name = "Linearized LPS bond-based"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage"
solver = "Quasi-Static"; % "Quasi-Static", "Dynamic/Explicit"
switch model.name
    case "PMB"
        alpha = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(horizon,omega);
        c1 = 6*E*1e6/mm;
        T = @interactionForce_PMB;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
    case "Linearized LPS"
        nu = 1/4;
        alpha = 1; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        mm = weightedVolume(horizon,omega); 
        lambda = E*nu/(1+nu)/(1-2*nu); mu = E/2/(1+nu);
        k = mu*(3*lambda + 2*mu)/(lambda + 2*mu);
        c1 = (2*k - 4*mu)*1e6/mm;
        c2 = 8*mu*1e6/mm;
        T = @interactionForce_LLPS;
        model.linearity = true;
        model.stiffnessAnal = false; % true if an analytical stiffness matrix for such model is implemented
        model.dilatation = true;
    case "Linearized LPS bond-based"
        alpha = 1; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        mm = weightedVolume(horizon,omega); 
        c1 = 6*E*1e6/mm;
        T = @interactionForce_LLPSBB;
        model.linearity = true;
        model.stiffnessAnal = true; % true if an analytical stiffness matrix for such model is implemented
        model.dilatation = false;
    case "Lipton Free Damage"
        %nu = 1/4;
        alpha = 1;
        mm = weightedVolume(horizon,omega);
        c1 = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2;
        c2 = (2*pi*horizon^3)^2/mm^2*E*1e6*(4*nu-1)/(2*(1+nu)*(1-2*nu))/2;
        T = @interactionForce_Lipton;
        model.linearity = true;
        model.stiffnessAnal = false;
        model.dilatation = true;
    case "LPS 2D"
        alpha = 1;
        mm = weightedVolume(horizon,omega);
        kappa = E*1e6/3/(1-2*nu); mu = E*1e6/2/(1+nu);
        c1 = kappa + mu/9*(nu+1)^2/(2*nu-1)^2;
        c2 = 8*mu/mm;
        T = @interactionForce_StateBased;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = true;
    otherwise
        error("Chosen model is not implemented or it was mistyped");
end

%% SIMULATION
tic
t_init = cputime;
for s_index = 1:length(sigma)
    % STRESS LOOP
    stresses = [0 sigma(s_index) 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for m_index = 1:1%length(m_vec)
        % HORIZON NUMBER LOOP
        %% -------------- Generate mesh -----------------
        h = h_vec(m_index); % grid spacing [m]
        m = m_vec(m_index); 
        a = 0.04/2; % height [m]
        b = 0.10/2; % length [m]
        [x,A] = generateMesh(h,[a b]); % Generates rectangular mesh 
        %% -------------- Boundary conditions ----------------
        [ndof,idb,bc_set,bodyForce,noFailZone] = boundaryCondition(x,stresses,m,h,A);
        %% -------------- GENERATE FAMILY ------------------
        [family,partialAreas,maxNeigh] = generateFamily_v2(x,horizon,m,m_index,true,"PA-AC"); % True for test
        %% -------------- Generate history variables ------------------
        history = historyDependency(x,maxNeigh);
        %% -------------- SOLVER -------------------
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t = 0:dt:40e-6; % 40 micro-secs simulation
                n_tot = length(t);
                [u_n,phi,energy] = solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,partialAreas,T,par_omega,history,noFailZone);
            case "Quasi-Static"
                n_tot = 2;
                [u_n,r,energy] = solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,T,par_omega,ndof,A);
            otherwise
                error("Solver not yet implemented.")
                pause
        end    
    end
end
toc
t_run = cputime - t_init;
%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        PostProcessing(x,u_n,n_tot,idb,energy,phi);
    case "Quasi-Static"
        PostProcessing(x,u_n,n_tot,idb,energy);
    otherwise
end