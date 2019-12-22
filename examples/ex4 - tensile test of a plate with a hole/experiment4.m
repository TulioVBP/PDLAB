function out = experiment4()
close all
clc
%% PARAMETERS
% --- Material --------
horizonVec = [0.005 0.004 0.003]; % [m]
lambda = 1.94e3; % [MPa]
G = 2.45e3; % [MPa]
E = G*(3*lambda + 2*G)/(lambda + G); % [MPa]
nu = lambda/2/(lambda+G);
rho = 2162;%3150*0.22 + 1442*0.66 + 1000*0.12; % [kg/m^3]
G0 = 2.28*10^3; % [J/m^2] [N/m]
for delta = 1:length(horizonVec)
horizon = horizonVec(delta);
% ---- Simulation ----------
sigmay = 6; % [MPa]
stress_app = '-';
m_vec = 4; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; alfa = 1;% Influence function options
par_omega = [horizon omega alfa];
notch_length = 10; % 1 cm
damage.crackIn = [-32.5 5;-32.5+notch_length 5]*10^-3; % Coordinates of the crack initial segment
b_parallelComp = false;
% ---- MODEL ---------
damage.damageOn = true; % True if applying damage to the model, false if not
model.name = "LPS 2D"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
solver = "Dynamic/Explicit"; % "Quasi-Static", "Dynamic/Explicit"
[model,c,T,damage] = models.modelParameters(model,par_omega,damage,E,nu,G0); % Check if it works    

%% SIMULATION
stresses = [0 sigmay 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
for m_index = 1:1%2:length(m_vec)
    % HORIZON NUMBER LOOP
    %% -------------- Generate mesh -----------------
    h = h_vec(m_index); % grid spacing [m]
    m = m_vec(m_index); 
    a = 0.120; % height [m]
    b = 0.065; % length [m]
    [x,A] = generateMesh_exp4(h,[a b]); % Generates rectangular mesh 
    %% -------------- Boundary conditions ----------------
    pc = prescribedBC(x,stresses); 
    [ndof,idb,bc_set,bodyForce,noFailZone] = mesh.boundaryCondition(x,stresses,m,h,A,3,stress_app,pc,damage);
    %% -------------- GENERATE FAMILY ------------------
    [family,partialAreas,maxNeigh,surfaceCorrection] = neighborhood.generateFamily_v2(x,A,horizon,m,m_index,false,"PA-HHB","None"); % True for test
    %% -------------- Generate history variables ------------------
    history = models.historyDependency(x,maxNeigh,model);
    %% -------------- SOLVER -------------------
    switch solver
        case "Dynamic/Explicit"
            dt = 0.5e-6; % 0.02 micro-sec is the one used by the paper
            t_tot = 1000e-6;% 10% of the actual simulation
            t = 0:dt:t_tot; % 40 micro-secs simulation
            n_tot = length(t);
            data_dump = 5;
            timeEval0 = tic;
            [t_s,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parallelComp,data_dump);
            t_cpu = toc(timeEval0);          
            filename = strcat('sim_m',int2str(m_vec(m_index)),'_d',int2str(horizon*1e3),'LPS','.mat');
            save(filename,'x','idb','u_n','phi','energy','t_cpu');
            out.x = x; out.idb = idb; out.un = u_n; out.energy = energy; out.t = t_cpu; out.phi = phi;
            disp(strcat('File saved as ',filename))
        case "Quasi-Static"
            if damage.damageOn
                error('Disable damage to run a quasi-static solver')
            end
            n_tot = 4;
            [u_n,r,energy] = solvers.solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,surfaceCorrection,T,c,model,par_omega,ndof,A,damage,history,noFailZone);
            out.x = x; out.un = u_n; out.energy = energy;
        otherwise
            error("Solver not yet implemented.")
            pause
    end
end

%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        PostProcessing_Exp4(x,u_n,n_tot,idb,energy,phi,t,t_s);
    case "Quasi-Static"
        postproc.PostProcessing(x,u_n,n_tot,idb,energy);
    otherwise
end
end
