function out = experiment4()
close all
clc
%% PARAMETERS
% --- Material --------
horizonVec = [0.005 0.004 0.002]; % [m]
lambda = 1.94e9; % [Pa]
G = 2.45e9; % [Pa]
E = G*(3*lambda + 2*G)/(lambda + G); % [MPa]
nu = lambda/2/(lambda+G);
rho = 2162;%3150*0.22 + 1442*0.66 + 1000*0.12; % [kg/m^3]
G0 = 2.28*10^3; % [J/m^2] [N/m]
for delta = 1:length(horizonVec)
    horizon = horizonVec(delta);
    % ---- Simulation ----------
    sigmay = 6e6; % [MPa]
    stress_app = '-';
    m_vec = 4; %[2 3 6 9]; % horizon number (last one should be 12 but we don't have enough space in memory) m = [2 3 6 12]
    h_vec = horizon./m_vec; % [m]
    omega = 3; alfa = 1;% Influence function options
    par_omega = [horizon omega alfa];
    notch_length = 10; % 1 cm
    damage.crackIn = [-32.5 5;-32.5+notch_length 5]*10^-3; % Coordinates of the crack initial segment
    b_parallelComp = false;
    % ---- MODEL ---------
    damage.damageOn = false; % True if applying damage to the model, false if not
    damage.DD = false;
    model.name = "LPS-T"; % "PMB", "Linearized LPS bond-based", "Lipton Free Damage" "LPS 2D" "Linearized LPS"
    solver = "Dynamic/Explicit"; % "Quasi-Static", "Dynamic/Explicit", "Quasi-Static Explicit"
    [modelo,damage] = models.modelParameters(model,par_omega,damage,E,nu,G0); % Check if it works    

    %% SIMULATION
    stresses = [0 sigmay 0]; % [sigma_x, sigma_y, tau_xy] - Pa/m^2
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
        values{1} = x; values{2} = maxNeigh; 
        modelo.history = values;  % Initialize his
        %% -------------- SOLVER -------------------
        switch solver
            case "Dynamic/Explicit"
                dt = 0.5e-6; % 0.02 micro-sec is the one used by the paper
                t_tot = 1000e-6;% 10% of the actual simulation
                t = 0:dt:t_tot; % 40 micro-secs simulation
                n_tot = length(t);
                data_dump = 5;
                timeEval0 = tic;
                [t_s,u_n,phi,energy] = solvers.solver_DynamicExplicit(x,t,idb,bodyForce,bc_set,family,A,partialAreas,surfaceCorrection,rho,modelo,par_omega,noFailZone,damage,b_parallelComp,data_dump);
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
                [u_n,r,energy] = solvers.solver_QuasiStatic(x,n_tot,idb,bodyForce,bc_set,family,partialAreas,surfaceCorrection,modelo,par_omega,ndof,A,damage,noFailZone);
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
                                                                               bc_set,family,A,partialAreas,surfaceCorrection,T,c,CC,model,par_omega,history,noFailZone,damage,b_parallelComp,beta,load_par,data_dump);
            otherwise
                error("Solver not yet implemented.")
                pause
        end
    end
end

%% POST-PROCESSING
switch solver
    case "Dynamic/Explicit"
        PostProcessing_Exp4(x,u_n,n_tot,idb,energy,phi,t,t_s);
    case "Quasi-Static"
        postproc.PostProcessing(x,u_n,n_tot,idb,energy);
    case "Quasi-Static Explicit"
        postproc.PostProcessing(x,u_n,index_s,idb,energy,phi);
    otherwise
end
end
