function out = experiment_template()
close all
clc
%% PARAMETERS
% --- Material --------
horizon = 0.05; % [m]
E = 72e3; % [MPa]
nu = 0.2;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
% --- PD ----------
m = 4; % mesh ration (m = horizon/h)
h = horizon./m; % [m]
omega = 3; gamma = 1;% Influence function options (1 - Exp., 2 - constant, 3 - conical)
par_omega = [horizon omega gamma];
PA_alg = "PA-AC"; % "FA", "PA-HBB", "PA-AC"
SE_alg = "None"; % "None", "Volume method"
% --- Mesh -----------------
a = 0.15; % height [m]
b = 1; % length [m]
[x,A] = mesh.generateMesh(h,[a b]); % Generates rectangular mesh 
% --- Initial damage ----
notch_length = 0.05; % Example 5 cm
damage.crackIn = [-0.3 -0.075;-0.3 -0.075+notch_length]; % Coordinates of the crack initial segment
% ---- MODEL ---------
damage.damageOn = false; % True if applying damage to the model, false if not
model.name = "PMB DTT"; % "PMB", "LBB", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
solver = "Quasi-Static"; % "Quasi-Static", "Dynamic/Explicit"
[model,c,T,damage] = models.modelParameters(model,par_omega,damage,E,nu,G0); % Check if it works    

%% SIMULATION
b_parallelComp = false; % true for parallel computation

% -------------- Boundary conditions ----------------
sigmay = 6; % [MPa] Example
astresses = [0 sigmay 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
stress_app = '-'; %'-' means BC provided directly (use prescribedBC.m). 't', 'l', 'r' and 'b' represents top, left, right and bottom edges of retangular domain (tbr in next updates)
pc = prescribedBC(x,stresses); % SET YOUR BCs IN THIS FUNCTION
[ndof,idb,bc_set,bodyForce,noFailZone] = mesh.boundaryCondition(x,stresses,m,h,A,3,stress_app,pc);

% -------------- GENERATE FAMILY ------------------
[family,partialAreas,maxNeigh,surfaceCorrection] = neighborhood.generateFamily_v2(x,A,horizon,m,1,false,PA_alg,SE_alg); % PA Algs: FA (choose for inhomogeneous mesh), PA-HHB and PA-AC. SE Algs: "None", "Volume"
        % It is possible to save multiple family files. Each file will be
        % identified as 'familyX.mat' where X is an integer. Here, X=1
% -------------- Generate history variables ------------------
history = models.historyDependency(x,maxNeigh,model); 

% -------------- SOLVER -------------------
switch solver
    case "Dynamic/Explicit"
        dt = 0.5e-6; % Time step
        t_tot = 500e-6;% Final time
        t = 0:dt:t_tot;
        n_tot = length(t);
        data_dump = 4; % Interval time steps to record outputs
        timeEval0 = tic;
        
        [t_s,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parallelComp,data_dump);
        
        t_cpu = toc(timeEval0);       
        filename = strcat('sim_m',int2str(m),'_d',int2str(horizon*1e3),'PMB','.mat'); % Choose your file name
        save(filename,'x','idb','u_n','phi','energy','t_cpu');
        out.x = x; out.idb = idb; out.un = u_n; out.energy = energy; out.t = t_cpu; out.phi = phi;
        disp("File saved as " +filename)
    case "Quasi-Static"
        if damage.damageOn
            error('Disable damage to run a quasi-static solver')
        end
        n_tot = 4; % Number of load steps
        [u_n,r,energy] = solvers.solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,surfaceCorrection,T,c,model,par_omega,ndof,A,damage,history,noFailZone);
        out.x = x; out.un = u_n; out.energy = energy;
    otherwise
        error("Solver not yet implemented.")
        pause
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
