% MAIN SCRIPT TO RUN A 2D PERIDYNAMIC BOND BASED MODEL - TEST
function out = experiment1_QS()
%% PARAMETERS
% --- Material --------
horizon = [5]*1e-3; % [m]
E = 72e9; % [MPa]
nu = 0.22;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
%damage.Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
% ---- Simulation ----------
sigmay = 4e6; % [MPa]
stress_app = 'tlr';
m_vec = 4; %[2 3 6 9]; % horizon number
h_vec = horizon./m_vec; % [m]
omega = 3; alfa = 1;% Influence function options
notch_length = 0.025; % 5 cm
%damage.crackIn = [0 0.01; notch_length 0.01]; % Coordinates of the crack initial segment
damage.crackIn = [];
b_parallelComp = false;
% ---- MODEL ---------
damage.damageOn = false; % True if applying damage to the model, false if not
damage.DD = false;
model.name = "LPS-T"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
solver = "Quasi-Static"; % "Quasi-Static", "Dynamic/Explicit"
%% SIMULATION
for s_index = 1:length(sigmay)
    % STRESS LOOP
    stresses = [0 sigmay(s_index) 0]; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for delta_index = 1:length(horizon)
        % HORIZON NUMBER LOOP
        par_omega = [horizon(delta_index) omega alfa];
        [model,c,T,damage] = models.modelParameters(model,par_omega,damage,E,nu,G0); % Check if it works    
        %% -------------- Generate mesh -----------------
        h = h_vec(delta_index); % grid spacing [m]
        m = m_vec(1); 
        a = 0.04/2; % height [m]
        b = 0.10/2; % length [m]
        [x,A] = mesh.generateMesh(h,[a b]); % Generates rectangular mesh 
        %% -------------- Boundary conditions ----------------
        [ndof,idb,bc_set,bodyForce,noFailZone] = mesh.boundaryCondition(x,stresses,m,h,A,true,stress_app);
        %% -------------- GENERATE FAMILY ------------------
        %[family,partialAreas,XJ,YJ,maxNeigh,surfaceCorrection] = generateFamily_v3(x,horizon,m,m_index,true,"None","None"); % True for test
        [family,partialAreas,maxNeigh,surfaceCorrection] = neighborhood.generateFamily_v2(x,A,horizon(delta_index),m,delta_index,false,"PA-AC","None"); % True for test
        %% -------------- Generate history variables ------------------
        history = models.historyDependency(x,maxNeigh,model);
        %% -------------- SOLVER -------------------
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t_tot = 30e-6;
                t = 0:dt:t_tot; % 40 micro-secs simulation
                n_tot = length(t);
                data_dump = 4; % Interval time steps to record outputs
                timeEval0 = tic;
                [t_s,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parallelComp,data_dump);
                t_cpu = toc(timeEval0);                      
                filename = strcat('sim_m',int2str(m_vec(1)),'_d',int2str(horizon(delta_index)*10^3),model.name,'.mat');
                save(filename,'x','idb','u_n','phi','energy','t_cpu');
                out.x = x; out.idb = idb; out.un = u_n; out.energy = energy; out.t = t_cpu;
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
end

%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        postproc.PostProcessing_Dyn(x,u_n,n_tot,idb,energy,phi,t,t_s);
    case "Quasi-Static"
        postproc.PostProcessing(x,u_n,n_tot,idb,energy);
    otherwise
end
end
