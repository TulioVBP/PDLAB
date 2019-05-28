% MAIN SCRIPT TO RUN A 2D PERIDYNAMIC BOND BASED MODEL - TEST
clear all
clc
close all

%% PARAMETERS
% --- Material --------
horizon = 4e-3; % [m]
E = 72e3; % [MPa]
nu = 1/3;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
damage.Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
% ---- Simulation ----------
sigma = 12; % [MPa]
m_vec = [2 4]; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; alfa = 1;% Influence function options
par_omega = [horizon omega alfa];
notch_length = 0.05; % 5 cm
damage.crackIn = [0 0.02; notch_length 0.02]; % Coordinates of the crack initial segment
b_parallelComp = false;
% ---- MODEL ---------
damage.damageOn = true; % True if applying damage to the model, false if not
model.name = "PMB"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
solver = "Dynamic/Explicit"; % "Quasi-Static", "Dynamic/Explicit"
switch model.name
    case "PMB"
        alfa = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(par_omega);
        c(1) = 6*E*1e6/mm;
        T = @interactionForce_PMB;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
        % Damage dependent Sc
        if true
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
    case "Linearized LPS bond-based"
        alfa = 0; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        mm = weightedVolume(par_omega); 
        c(1) = 6*E*1e6/mm;
        T = @interactionForce_LLPSBB;
        model.linearity = true;
        model.stiffnessAnal = true; % true if an analytical stiffness matrix for such model is implemented
        model.dilatation = false;
    case "Lipton Free Damage"
        alfa = 1;
        mm = weightedVolume(par_omega);
        c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2;
        c(2) = (2*pi*horizon^3)^2/mm^2*E*1e6*(4*nu-1)/(2*(1+nu)*(1-2*nu))/2;
        T = @interactionForce_Lipton;
        model.linearity = true;
        model.stiffnessAnal = false;
        model.dilatation = true;
    case "LPS 2D"
        alfa = 1;
        mm = weightedVolume(par_omega);
        kappa = E*1e6/3/(1-2*nu); mu = E*1e6/2/(1+nu);
        c(1) = kappa + mu/9*(nu+1)^2/(2*nu-1)^2;
        c(2) = 8*mu/mm;
        c(3) = nu;
        T = @interactionForce_StateBased;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = true;
%     case "Linearized LPS"
%         %nu = 1/4;
%         alfa = 1; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
%         mm = weightedVolume(par_omega); 
%         lambda = E*nu/(1+nu)/(1-2*nu); mu = E/2/(1+nu);
%         k = mu*(3*lambda + 2*mu)/(lambda + 2*mu);
%         c(1) = (2*k - 4*mu)*1e6/mm;
%         c(2) = 8*mu*1e6/mm;
%         T = @interactionForce_LLPS;
%         model.linearity = true;
%         model.stiffnessAnal = false; % true if an analytical stiffness matrix for such model is implemented
%         model.dilatation = true;
    otherwise
        error("Chosen model is not implemented or it was mistyped");
end

%% SIMULATION
for s_index = 1:length(sigma)
    % STRESS LOOP
    stresses = [0 sigma(s_index) 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for m_index = 1:1%length(m_vec)
        % HORIZON NUMBER LOOP
        %% -------------- Generate mesh -----------------
        h = h_vec(m_index); % grid spacing [m]
        m = m_vec(m_index); 
        a = 0.04; % height [m]
        b = 0.10; % length [m]
        [x,A] = generateMesh(h,[a b]); % Generates rectangular mesh 
        %% -------------- Boundary conditions ----------------
        [ndof,idb,bc_set,bodyForce,noFailZone] = boundaryCondition(x,stresses,m,h,A,true,'tlr');
        %% -------------- GENERATE FAMILY ------------------
        [family,partialAreas,maxNeigh,surfaceCorrection] = generateFamily_v2(x,horizon,m,m_index,true,"PA-AC","None"); % True for test
        %% -------------- Generate history variables ------------------
        history = historyDependency(x,maxNeigh,model);
        %% -------------- SOLVER -------------------
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t_tot = 12e-6;%40e-6;
                t = 0:dt:t_tot; % 40 micro-secs simulation
                n_tot = length(t);
                tic
                [u_n,phi,energy] = solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parallelComp);
                toc
            case "Quasi-Static"
                if damage.damageOn
                    error('Disable damage to run a quasi-static solver')
                end
                n_tot = 2;
                [u_n,r,energy] = solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,surfaceCorrection,T,c,model,par_omega,ndof,A,damage);
            otherwise
                error("Solver not yet implemented.")
                pause
        end    
    end
end

%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        PostProcessing(x,u_n,n_tot,idb,energy,phi);
    case "Quasi-Static"
        PostProcessing(x,u_n,n_tot,idb,energy);
    otherwise
end
