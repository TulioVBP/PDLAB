function [t_s,u_load,un_sample,index_s,phi_sample,energy,history,time_up,F_load,CC,damage] = solver_QuasiStaticExplicit(x,idb,body_force,...
   bc_set,familyMat,A,partialAreas,surfaceCorrection,damping,model,par_omega,noFailZone,damage,b_parll,beta,load_par,data_dump,b_partialSim,perc_ps)
% Unravel load parameters
n_load = load_par.n_load;
n_iterMax = load_par.n_iterMax;

% Explicit time solver with dynamic relaxation for quasi-static experiments
   if nargin < 17 % No data dump
        data_dump = 1;
        b_partialSim = [false false];
   elseif nargin < 18
        b_partialSim = [false false];
   end
   if nargin < 19
       perc_ps = 1;
   end
    %% Create volume
    h = norm(x(1,:)-x(2,:));
    if length(A)== 1
        A = A*ones(length(x),1); % homogeneous volume
    end
    % {No fail to damage variable}
    damage.noFail = noFailZone;   
    
    %% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
    
    %% Defining cracking trespassing matrix
    crackSegments = size(damage.crackIn,1); % At least 2
    damage.checkCrack = zeros(size(model.history.S));
    for ii = 1:size(model.history.S)
        x_j = x(familyMat(ii,familyMat(ii,:)~=0),:);
        check = zeros(size(x_j,1),crackSegments-1);
        for kk = 1:crackSegments-1
            for jj = 1:size(x_j,1)
            [~,check(jj,kk)] = neighborhood.checkBondCrack(x(ii,:),x_j(jj,:),damage.crackIn(kk,:),damage.crackIn(kk+1,:));
            end
        end
        check = check';
        if isempty(check)
           brokenBonds = logical(zeros(size(x_j,1),1));
        elseif size(check,1) == 1
           brokenBonds = logical(check);
        else
           brokenBonds = any(check);
        end
        damage.brokenBonds(ii,1:length(x_j)) = brokenBonds; 
    end
    disp('Check for broken bonds done.')
    
    %% INITIALIZE SIMULATION MATRICES
    ndof = 2*length(x) - size(bc_set,1); % Number of degree of freedoms
    u_load = zeros(2*length(x),n_load); % Displacement for load steps
    mu = cell(length(x),n_load); % Damage factor initialization
    damage.phi = zeros(length(x),1); % Initializing phi in the damage structure
    history.S = model.history.S; %
    if model.b_dilatation
        history.theta = model.history.theta;
    end
    
    % Initialize time vectors and density according to this
    if load_par.b_tgiven
        t = load_par.t;
        dt = abs(t(2)-t(1));
        % Density
        u0 = zeros(2*length(x),1);
        lambda = evaluateLambda(x,u0,ndof,idb,familyMat,partialAreas,surfaceCorrection,A,par_omega,model,damage,history,dt,mu(:,1));
        lambda_inv = 1./lambda;
        figure
        histogram(lambda)
    else
        % Density
        rho = load_par.rho;
        if length(rho) == 1
            lambda = rho*ones(ndof,1);
        else
            lambda_full = zeros(length(rho),1);
            for ii = 1:length(rho)
                lambda_full(dof_vec(ii,:)) = rho(ii);
            end
            lambda = lambda_full(1:ndof);
        end
        lambda_inv = 1./lambda;
        
        % Time vector
        u0 = zeros(2*length(x),1);
        dt = evaluateDeltaT(x,u0,ndof,idb,familyMat,partialAreas,surfaceCorrection,A,par_omega,model,damage,history,lambda,mu(:,1));
        t = 0:dt:load_par.t_max;
    end
    t_full = t(1):dt:t(end)*n_load;
    t_full = [t_full t_full+dt:dt:(t_full+n_iterMax*dt)];
    
    % Initial condition
    u_n = zeros(2*length(x),2); % 2*N x 2  2D matrix (length(t_full))
    v_n = zeros(2*length(x),2); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.
    phi = zeros(length(x),length(t_full)); % No initial damage
    
    % Checking boundary conditions
    if size(body_force,2) ~= 1
        error('Prescribed body forces must be constant.')
    end
    
    if sum(bc_set(:,3) ~= 0 & bc_set(:,2) ~= 0)
        error('The boundary condition matrix bc_set has prescribed both displacement and velocity for the same node.')
    end
                
    if sum(bc_set(:,3)) ~= 0
        warning('No velocity constraints should be added for this solver. The solver will ignore them.')
        bc_set = bc_set(bc_set(:,3)==0,:); % Eliminating the dofs with velocity constraint
    end
    
    
    % {Defining n_sample and initial energy variables
    n_sample = [1 data_dump:data_dump:length(t_full)];
    if n_sample(end)~=length(t_full)
        n_sample = [n_sample length(t_full)];
    end
    t_s = t_full(n_sample); % Sample time
    % Maybe the above is not necessary
    index_s = 0;
    
    energy.W = zeros(length(x),length(n_sample)); % No initial  deformation hence no initial strain energy density
    energy.KE = zeros(length(x),length(n_sample));
    energy.EW = zeros(length(x),length(n_sample));
    un_sample = zeros(2*length(x),length(n_sample));
    energy_ext_var = zeros(length(x),1); % External energy that takes into account variable body force (velocity constraints)
    fn = zeros(2*length(x),1); % Initial force
    u_const = zeros(length(v_n)-(ndof),1); % Constraint nodes
    theta = zeros(length(x),1); % Preallocate theta
    
    F_load = zeros(length(n_sample),2);
    phi_sample = zeros(length(x),length(n_sample));
    
    % Temporary variables
    history_tempS = history.S;
    if model.b_dilatation 
        history_tempT = history.theta;
    end
    fn_temp = zeros(size(x));
    phi_temp = zeros(length(x),1);
    
    n_initial = 1;
    kk = 0;
    bc_set_part = bc_set;
    crit_var = zeros(1,length(t)-1);
    c_Kneg = 0;
    C = 0;
    
    % {Recoverying temporary files}
    n_final_load = floor(n_load*perc_ps);
    if exist('tempsim.mat','file') && b_partialSim(1)
        load('tempsim.mat');
        n_initial = n;
        disp('Found a partial simulation. Continuing it...')
    end
    
    timerVal = tic; % Initiating another tic
    b_crack = true;
    b_load = true;
    b_finalIter = false;
    
    
    while b_load
        %% #########  LOAD LOOP ##########
        kk = kk+1*(~b_finalIter);
        bc_set_part(:,2) = bc_set(:,2)*kk/n_load; % Upload boundary condition
        body_force_part(:,1) = body_force*kk/n_load; % bn
        % Obtain data from previous simulation
        if kk >1
            mu(:,kk) = mu(:,kk-1);
            n_initial = n+1;
        end
        F = zeros(2*length(x),2);
        
        n_in = n_initial; n_fin = (n_initial+length(t)*(~b_finalIter)+n_iterMax*(b_finalIter)-2*(~b_finalIter));
        
        for n = n_in:n_fin
            % Instatiate body force
            bn = body_force_part; % Increment b
           
            %% ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % #### Step 0 - Acceleration at a_n
            an =  lambda_inv .* (fn(1:ndof) + bn(1:ndof) - C * lambda.* v_n(1:ndof,1));
            if sum(isnan(an)) || sum(isnan(fn(1:ndof))) || sum(sum(isnan(v_n(1:ndof,:))))
                disp('an is nan')
            end
            % #### Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt/2 * an; % V(n+1/2)
            % #### Step 2 - Update displacement
            u_n(:,1) = u_n(:,2); % u(n) is stored in the u(1)
            u_n(1:ndof,2) = u_n(1:ndof,1) + dt*v_n(1:ndof,2); % u(n+1) - %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
             % ----- Solving for the displacement constraint nodes ----
            if ~isempty(u_const)
                u_const = bc_set_part(:,2); % Defining the displacements for the nodes with no velocity
                u_n(ndof+1:end,2) = u_const;
            end
            u1 = u_n(:,1);
            u2 = u_n(:,2); % Vector of displacement that will be used to come back to 2D

            % ---- {Evaluating dilatation} ----
           
            damage.phi = phi(:,n); % Accessing current damage situation
            if model.b_dilatation
                if b_parll
                    parfor ii = 1:length(x)
                        [theta(ii),history_tempT(ii)] = model.dilatation(x,u2,familyMat(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,damage,history.S, history.theta);
                    end
                else
                    [theta,history_tempT] = model.dilatation(x,u2,familyMat,partialAreas,surfaceCorrection,[],idb,par_omega,damage,history.S,history.theta);
                end
                history.theta = history_tempT; % Assigning up-to-date history variable
            end
            
            % ####### Step 3 - Update velocity 
            % Storing fn
            F(:,1) = fn + bn; % F(n)
            
            % {Evaluate f[n+1]}
            energy_pot = zeros(length(x),1); % Pre-allocate potential energy
            b_sampling = rem(n+1,data_dump) == 0 || (n+1 == n_fin && b_finalIter); % Deciding when to evaluate the energy

            if b_parll
                parfor ii = 1:length(x)
                   [fn_temp(ii,:),history_tempS(ii,:),mu{ii,kk},phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,dt,history,A,theta,b_sampling);
                end
            else 
                for ii = 1:length(x)
                   [fn_temp(ii,:),history_tempS(ii,:),mu{ii,kk},phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,dt,history,A,theta,b_sampling);

                end
            end
                
            % Converting the temporary variables
            for ii = 1:length(x)
                fn(dof_vec(ii,:)) = fn_temp(ii,:)';
            end
            history.S = history_tempS; % Updating the history variable related to the stretch
            phi(:,n+1) = phi_temp;
      
            % Evaluate V(n+1)
            an =  lambda_inv .* (fn(1:ndof) + bn(1:ndof) - C * lambda.* v_n(1:ndof,2)); % A(n+1)
            v_n(1:ndof,1) = v_n(1:ndof,2) + dt/2*an; % V(n+1) is stored in the next V(n)
            
            % ####### Step 5 - Determine C
            F(:,2) = fn + bn; % F(n + 1)
            if damping == 'critical'
                [C,c_Kneg] = evaluateDamping(lambda,u_n(:,2),ndof,dt,F,v_n(:,2),c_Kneg);
            else
                C = damping;
            end
            CC(n) = C;      
            
            %% Sampling quantities
            if b_sampling
                index_s = index_s + 1;%n_sample == n+1;
                % {Un sampled}
                un_sample(:,index_s) = u_n(:,2);
                
                % {Force load} - EDIT HERE TO OBTAIN FORCES CROSSING A GIVEN TRANSVERSAL SECTION 
                F_load(index_s,:) = forceDisplacement(x,u_n(:,2),familyMat,dof_vec,partialAreas,surfaceCorrection,par_omega,model,damage,history,theta,A);
                
                % {Phi}
                phi_sample(:,index_s) = phi_temp;
                
                % {Potential energy}
                energy.W(:,index_s) = energy_pot;
                
                % {External work}
                % Constant body force
                energy_ext = dot(u2(dof_vec)',bn(dof_vec)')'.*A;
                
                % {External work realized by the velocity constraint}
                if ~isempty(bc_set)
                    vel_dof = idb(bc_set(bc_set(:,2)~=0,1));
                    u_e1 = u_n(vel_dof,1);
                    u_e2 = u_n(vel_dof,2);
                    du = (u_e2-u_e1);
                else
                    vel_dof = [];
                    du = [];
                end
                bf = -( F(vel_dof,2) + F(vel_dof,1) - 2*bn(vel_dof) )/2;
                add_ext = bf.*du;
                for gg = 1:length(vel_dof)
                    ind_vel = find(dof_vec == vel_dof(gg));
                    ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
                    energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
                end
                energy.EW(:,index_s) = energy_ext + energy_ext_var;
                
                % {Kinectic energy}
                for ll = 1:length(x)
                    dofk = dof_vec(ll,:);
                    if dofk(1) <= ndof
                        rho(1) = lambda(dofk(1));
                    else
                        rho(1) = 0;
                    end
                   if dofk(2) <= ndof
                        rho(2) = lambda(dofk(2));
                    else
                        rho = 0;
                    end
                    energy.KE(ll,index_s) =  sum(1/2*rho*v_n(dofk,1).^2.*A(ll));
                end
            else
                 % {External incremental work only}
                
                du = u2(dof_vec) - u1(dof_vec);
                energy_ext_var = energy_ext_var+dot(du',bn(dof_vec)')'.*A;

                % {External work realized by the velocity constraint}
                if ~isempty(bc_set)
                    vel_dof = idb(bc_set(bc_set(:,2)~=0,1));
                    u_e1 = u_n(vel_dof,1);
                    u_e2 = u_n(vel_dof,2);
                    du = (u_e2-u_e1);
                else
                    vel_dof = [];
                    du = [];
                end
                bf =-( F(vel_dof,2) + F(vel_dof,1) - 2*bn(vel_dof) )/2;
                add_ext = bf.*du;
                for gg = 1:length(vel_dof)
                    ind_vel = find(dof_vec == vel_dof(gg));
                    ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
                    energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
                end
            end

            %% ############ Verifying convergence ##########
            PHI = fn(1:ndof) + bn(1:ndof);
            if norm(bn(1:ndof)) ~= 0
                crit_var(n) = norm(PHI)/norm(bn(1:ndof));
            else
                crit_var(n) = 1;
            end
            if crit_var(n) < beta
                    disp("Convergence achieved for the load step " + num2str(kk) + " ...")
                    disp("Number of negative K's = " + int2str(c_Kneg))
                    if kk == n_final_load
                        b_load = false;
                    end
                    break
            elseif sum(phi(:,n+1) > phi(:,n)) && n > 2 && b_crack
                disp("Crack nucleation (propagation) observed at time step " + int2str(n))
                b_crack = false;
                disp('Ignore break and propagate simulation.')
            else
                 %% ############ COUNTING THE PROCESSING TIME #############
                time_up = toc(timerVal);
                if n >1
                   clc
                end
                %lineLength = fprintf("Load step " + num2str(kk) + " out of " + num2str(n_load) +": Time step "+ num2str(n) + "%. ETA: "+ num2str(time_up/n*(length(t)-n)));
                disp("Load step " + num2str(kk) + " out of " + num2str(n_load) +": Time step "+ num2str(n) + " out of " + num2str(length(t_full)) + ". Criterion: " + crit_var(n) +...
                    ". ETA: "+ num2str(time_up/n*(length(t_full)+n_iterMax-n)) + " (avg. per timestep : " + num2str(time_up/n) + "). C = " + num2str(C))        
            end
            %% Conditions to finish simulation if load is applied through disp.
            
            if n == n_fin && b_finalIter
                b_load = false; % Stop simulation
            end
            
            %% EDIT HERE YOUR CONDITION FOR FINAL ITERATION
            if n == n_fin && ((~b_crack && abs(F_load(index_s,2)) < max(F_load(:,2))*1/6) || kk == n_final_load) 
                b_finalIter = true;
            end
            
            if ~damage.damageOn && kk == n_final_load && n == n_fin
                 b_load = false;
            end
            
            %% Verifying rate of convergence
            if n > 2
                if (abs(crit_var(n) - crit_var(n-1)))/crit_var(n-1) < 10^-10 && norm(bn(1:ndof)) ~= 0
                    error("Convergence rate is too slow and convergence criteria was not met. Conv. crit. = " + num2str(crit_var(n)));
                end
            end
        end
        u_load(:,kk) = u_n(:,1);
    end
    if b_partialSim(2)
        clear n_final_load;
        save('tempsim.mat')
    end
    [phi_sample,energy,un_sample,F_load] = cap_sample(phi_sample,energy,un_sample,F_load,index_s); % Eliminating last zeros, if any
end

%%
function [f_i,history_upS,mu_j,phi_up,energy_pot] = parFor_loop(x,u_n,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,dt,history,A,theta,b_sampling)
   % Loop on the nodes
   family = familyMat(ii,familyMat(ii,:)~=0);
   history_upS = history.S(ii,:);
   neig_index = 1:length(family);
   jj = family(neig_index);
   
   % Loop on their neighbourhood
   noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
   Vj = partialAreas(ii,neig_index)'; 
   if model.b_dilatation
      [fij,history_upS(neig_index),mu_j] = model.T(x,u_n,theta,ii,jj,dof_vec,par_omega,Vj,damage,history.S(ii,neig_index),history.theta,noFail);
   else
      [fij,history_upS(neig_index),mu_j] = model.T(x,u_n,ii,jj,dof_vec,par_omega,Vj,damage,history.S(ii,neig_index),noFail);
   end
   lambda = surfaceCorrection(ii,neig_index)';
   f_i = sum(fij.*Vj.*lambda);
   % Damage index
   areaTot = sum(Vj);
   partialDamage = sum(mu_j.*Vj);
   phi_up = 1 - partialDamage/areaTot;
   if b_sampling
       if ~model.b_dilatation
           % Strain energy
           W = model.strainEnergyDensity(x,u_n,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,damage,history_upS(1,neig_index));
       else
           W = model.strainEnergyDensity(x,u_n,theta,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,damage,history_upS(neig_index),history.theta); % neig_index == length(family)
       end
       % Stored strain energy
       energy_pot = W.*A(ii);
   else
       energy_pot = 0;
   end
end

%% Reducing the components
function [xs] = sampling(x,t,ts)
    xs = zeros(size(x,1),length(ts));
    for iii = 1:size(x,1)
        xs(iii,:) = interp1(t,x(iii,:),ts);
    end
end

function lambda = evaluateLambda(x,u,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,model,damage,history,dt,mu)
    % 0 - Correct mu
    if isempty(mu{1})
        mu(:) = {1};
    end
    % 1 - Define stiffness matrix
    if model.b_stiffnessAnal
        K = model.analyticalStiffnessMatrix(x,u,ndof,idb,family,partialAreas,surfaceCorrection,ones(size(V)),par_omega,damage,history,mu);
    else
        K = tangentStiffnessMatrix(x,u,idb,family,partialAreas,ones(size(V)),surfaceCorrection,ndof,par_omega,model,damage,history);
    end
    % 2 - Define lambda
    lambda = diag(eye(ndof));
    for ii = 1:ndof
        lambda(ii) = lambda(ii) * (2/4*dt^2*sum(abs(K(ii,1:ndof))));
    end
    lambda(lambda < 1e-13) = max(lambda(lambda > 1e-13));
end

function dt = evaluateDeltaT(x,u,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,model,damage,history,lambda,mu)
    % 0 - Correct mu
    if isempty(mu{1})
        mu(:) = {1};
    end
    % 1 - Define stiffness matrix
    if model.b_stiffnessAnal
        K = model.analyticalStiffnessMatrix(x,u,ndof,idb,family,partialAreas,surfaceCorrection,ones(size(V)),par_omega,damage,history,mu);
    else
        K = tangentStiffnessMatrix(x,u,idb,family,partialAreas,ones(size(V)),surfaceCorrection,ndof,par_omega,damage,history);
    end
    % 2 - Define lambda
    dt_vec = diag(eye(ndof));
    for ii = 1:ndof
        dt_vec(ii) = sqrt(4*lambda(ii)/2/sum(abs(K(ii,1:ndof))));
    end
    dt = min(dt_vec);
    figure
    histogram(dt_vec,30);
end

function [C,count] = evaluateDamping(lambda,u,ndof,dt,F,v,count)
    % 4 - Define K diagonal
    vv = v(1:ndof);
    Kn = -(F(1:ndof,2) - F(1:ndof,1))./lambda./(dt*vv);
    Kn (vv < 1e-14 & vv > -1e-14) = 0;
    Kn = diag(Kn);
    % 5 - Evaluate the damping
    if u(1:ndof)'*Kn*u(1:ndof) < 0
        count = count + 1;
    end
    C = abs(2*sqrt(u(1:ndof)'*Kn*u(1:ndof)/(u(1:ndof)'*u(1:ndof))));
    disp("C is "+num2str(C));
    if C < 0 || isnan(C)
        C = 0;
    end
end

function F = forceDisplacement(x,u,family,dof_vec,partialAreas,surfaceCorrection,par_omega,model,damage,history,theta,V)
   x_bound = 0.020; 
   F = [0 0];
   for ii = 1:length(x)
       if x(ii,2) < x_bound-1e-13
              family_ii_comp = family(ii,family(ii,:)>0);
              neig_index = 1:length(family_ii_comp);
              family_ii = family_ii_comp(x(family_ii_comp,2) > x_bound);
              partialAreas_ii = partialAreas(ii,x(family_ii_comp,2) > x_bound);
              surfaceCorrection_ii = surfaceCorrection(ii,x(family_ii_comp,2) > x_bound);
              neig_index = neig_index(x(family_ii_comp,2) > x_bound);
              
              if ~isempty(family_ii)
                f_i = forceSection(x,u,dof_vec,ii,family_ii,neig_index, partialAreas_ii,surfaceCorrection_ii,par_omega,model,damage,history,theta);
                F = F + f_i*V(ii);
              end
       end
   end
end

function f_i = forceSection(x,u_n,dof_vec,ii, family, neig_index,partialAreas,surfaceCorrection,par_omega,model,damage,history,theta)
   jj = family();
   % Loop on their neighbourhood
   noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
   Vj = partialAreas';
   if model.b_dilatation
      [fij,~,~] = model.T(x,u_n,theta,ii,jj,dof_vec,par_omega,Vj,damage,history.S(ii,neig_index),history.theta,noFail);
   else
      [fij,~,~] = model.T(x,u_n,ii,jj,dof_vec,par_omega,Vj,damage,history.S(ii,neig_index),noFail);
   end
   lambda = surfaceCorrection';
   
   f_i = sum(fij.*Vj.*lambda);
   
end

function [phi_sample,energy,un_sample,F_load] = cap_sample(phi_sample,energy,un_sample,F_load,index_s)
    phi_sample = phi_sample(:,1:index_s);
    energy.EW = energy.EW(:,1:index_s);
    energy.KE = energy.KE(:,1:index_s);
    energy.W = energy.W(:,1:index_s);
    un_sample = un_sample(:,1:index_s);
    F_load = F_load(1:index_s,:);
end