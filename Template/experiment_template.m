function out = experiment_template()
close all
clc
%% PARAMETERS
% --- Material --------
horizon = 0.2; % [m]
E = 72e9; % [Pa]
nu = 0.2;
rho = 2440; % [kg/m^3]
G0 = 3.8; % [J/m^2]
% --- PD ----------
m = 4; % mesh ration (m = horizon/h)
h = horizon./m; % [m]
omega = 3; gamma = 1;% Influence function options (1 - Exp., 2 - constant, 3 - conical)
par_omega = [horizon omega gamma];
PA_alg = "PA-HHB"; % "FA", "PA-HHB", "PA-AC"
SE_alg = "None"; % "None", "Volume method"
dt = 1;%0.5e-6; % Time step

% --- Mesh -----------------
a = 0.15; % height [m]
b = 1; % length [m]
[x,A] = mesh.generateMesh(h,[a b]); % Generates rectangular mesh 

% --- Initial damage ----
notch_length = 0.05; % Example 5 cm
damage.crackIn = [-0.3 -0.075;-0.3 -0.075+notch_length]; % Coordinates of the crack initial segment
damage.DD = false; % Damage dependent criteria

% ---- MODEL ---------
damage.damageOn = true; % True if applying damage to the model, false if not
model.name = "LSJ-T"; % "PMB", "DTT", "LBB", "LSJ-T", "LPS-T", "Linearized LPS"
solver = "Quasi-Static Explicit"; % "Quasi-Static", "Dynamic/Explicit", "Quasi-Static Explicit"
[model,damage,modelo] = models.modelParameters(model,par_omega,damage,E,nu,G0,dt); % Check if it works    

%% SIMULATION
b_parallelComp = false; % true for parallel computation

% -------------- Boundary conditions ----------------
sigmay = 6; % [MPa] Example
stresses = [0 sigmay 0]*1e6; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
stress_app = '-'; %'-' means BC provided directly (use prescribedBC.m). 't', 'l', 'r' and 'b' represents top, left, right and bottom edges of retangular domain (tbr in next updates)
pc = prescribedBC(x,stresses); % SET YOUR BCs IN THIS FUNCTION
[ndof,idb,bc_set,bodyForce,noFailZone] = mesh.boundaryCondition(x,stresses,m,h,A,3,stress_app,pc,damage);

% -------------- GENERATE FAMILY ------------------
[family,partialAreas,maxNeigh,surfaceCorrection] = neighborhood.generateFamily_v2(x,A,horizon,m,1,false,PA_alg,SE_alg); % PA Algs: FA (choose for inhomogeneous mesh), PA-HHB and PA-AC. SE Algs: "None", "Volume"
        % It is possible to save multiple family files. Each file will be
        % identified as 'familyX.mat' where X is an integer. Here, X=1
% -------------- Generate history variables ------------------
values{1} = x; values{2} = maxNeigh; 
modelo.history = values;  % Initialize his

% -------------- SOLVER -------------------
switch solver
    case "Dynamic/Explicit"
        t_tot = 1000e-6;% Final time
        t = 0:dt:t_tot;
        n_tot = length(t);
        data_dump = 4; % Interval time steps to record outputs
        timeEval0 = tic;
        
        [t_s,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,rho,modelo,par_omega,noFailZone,damage,b_parallelComp,data_dump);
        
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
    case "Quasi-Static Explicit"
        load_par.n_iterMax = 50;
        load_par.n_load = 100; % Number of load steps
        CC = 'critical'; % Change to a numerical value if you want CC fixed
        data_dump = 8; % Interval time steps to record outputs
        beta = 5e-6; % Useful for convergence criteria
        if true % CHANGE IF YOU WANT TO INPUT THE DENSITY INSTEAD
            % TIME GIVEN
            load_par.b_tgiven = true;
            dt = 1; % Time step
            t_tot = 30;% Final time
            load_par.t = 0:dt:t_tot;
        else
            % DENSITY GIVEN
            lambda = rho; %Alternatively, you can save a previous virtual density matrix and then load('lambda.mat','lambda'). It is usefull in convergence studies where the mesh changes.
            load_par.b_tgiven = false;
            load_par.rho = max(lambda);
            load_par.t_max = 30;
        end
        [t_s,u_load,u_n,index_s,phi,energy,history,time_up,F_load,CC] = solvers.solver_QuasiStaticExplicit(x,idb,bodyForce,...
                                                                       bc_set,family,A,partialAreas,surfaceCorrection,CC,modelo,par_omega,noFailZone,damage,b_parallelComp,beta,load_par,data_dump);
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
    case "Quasi-Static Explicit"
        postproc.PostProcessing(x,u_n,index_s,idb,energy,phi);
    otherwise
end
end
