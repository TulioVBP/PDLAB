function out = experiment2_Dyn()
close all
clc
%% PARAMETERS
% --- Material --------
horizon = 6e-3; % [m]
E = 72e9; % [Pa]
nu = 0.22;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
%damage.Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon); % Critical elongation - conical micromodulus
% ---- Simulation ----------
sigmay = 2e6; % [Pa]
stress_app = 'tblr';
m_vec = [2 4]; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
h_vec = horizon./m_vec; % [m]
omega = 3; alfa = 1;% Influence function options
par_omega = [horizon omega alfa];
notch_length = 0.05; % 5 cm
damage.crackIn = [-0.05 0; notch_length-0.05 0]; % Coordinates of the crack initial segment
b_parallelComp = false;
% ---- MODEL ---------
damage.damageOn = true; % True if applying damage to the model, false if not
damage.DD = false;
model.name = "LSJ-T"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
solver = "Dynamic/Explicit"; % "Quasi-Static", "Dynamic/Explicit"
[model,c,T,damage] = models.modelParameters(model,par_omega,damage,E,nu,G0); % Check if it works    

%% SIMULATION
for s_index = 1:length(sigmay)
    % STRESS LOOP
    stresses = [0 sigmay(s_index) 0]; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
    for m_index = 2:2%2:length(m_vec)
        % HORIZON NUMBER LOOP
        %% -------------- Generate mesh -----------------
        h = h_vec(m_index); % grid spacing [m]
        m = m_vec(m_index); 
        a = 0.04; % height [m]
        b = 0.10; % length [m]
        [x,A] = mesh.generateMesh(h,[a b]); % Generates rectangular mesh 
        %% -------------- Boundary conditions ----------------
        [ndof,idb,bc_set,bodyForce,noFailZone] = mesh.boundaryCondition(x,stresses,m,h,A,false,stress_app,[],damage);
        %% -------------- GENERATE FAMILY ------------------
        [family,partialAreas,maxNeigh,surfaceCorrection] = neighborhood.generateFamily_v2(x,A,horizon,m,m_index,false,"PA-AC","None"); % True for test
        %% -------------- Generate history variables ------------------
        history = models.historyDependency(x,maxNeigh,model);
        %% -------------- SOLVER -------------------
        switch solver
            case "Dynamic/Explicit"
                dt = 0.02e-6; % 0.02 micro-sec is the one used by the paper
                t_tot = 43e-6;% 10% of the actual simulation
                t = 0:dt:t_tot; % 40 micro-secs simulation
                n_tot = length(t);
                data_dump = 4; % Interval time steps to record outputs
                timeEval0 = tic;
                [ts,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parallelComp,data_dump);
                t_cpu = toc(timeEval0);       
                %filename = strcat('sim_m',int2str(m_vec(m_index)),'_d',num2str(horizon),int2str(IDN(1)),'.mat');
                %filename = tempname;
                %IDN = randi([1 100],1,3);
                IDN = t_cpu;
                filename = strcat('sim_m',int2str(m_vec(m_index)),'_d',int2str(horizon*1e3),'PMB','.mat');
                save(filename,'x','idb','u_n','phi','energy','t_cpu');
                out.x = x; out.idb = idb; out.un = u_n; out.energy = energy; out.t = t_cpu;
                disp(strcat('File saved as ',filename))
            case "Quasi-Static"
                if damage.damageOn
                    error('Disable damage to run a quasi-static solver')
                end
                n_tot = 2;
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
