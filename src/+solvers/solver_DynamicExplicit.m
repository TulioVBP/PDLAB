% Explicit time - dynamic solver

function [t_s,u_n,phi,energy,history,time_up] = solver_DynamicExplicit(x,t,idb,body_force,bc_set,familyMat,A,partialAreas,surfaceCorrection,rho,model,par_omega,noFailZone,damage,b_parll,data_dump)
    if nargin < 19 % No data dump
        data_dump = 1;
    end
    ndof = 2*length(x) - size(bc_set,1);
    %% Create volume
    h = norm(x(1,:)-x(2,:));
    if length(A)== 1
        A = A*ones(length(x),1); % homogeneous volume
    end
    % {No fail to damage variable}
    damage.noFail = noFailZone;   
    %% Verify if dt is small enough
     dt = abs(t(2) - t(1));
%     dt_crit = criticalTimeStep(x,familyMat,partialAreas,par_omega,c,rho,model);
%     dt_ratio = dt/dt_crit;
%     if dt_ratio < 1 && dt_crit > 0
%         disp("Time-step " + num2str(dt) +" sec < critical time-step, " + num2str(dt_crit) + ...
%             " sec. Safety factor: " + num2str(dt_ratio) + ". The simulation should converge.")
%     elseif dt > 0
%        disp("Time-step " + num2str(dt) + " sec > critical time-step, " + num2str(dt_crit) + ...
%             " sec. Safety factor: " + num2str(dt_ratio)+ ". The simulation shall explode.")
%     end
    %% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
    %% {Defining cracking trespassing matrix}
    crackSegments = size(damage.crackIn,1); % At least 2
    damage.checkCrack = zeros(size(model.history.S));
    for ii = 1:size(model.history.S,1)
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
        Minv = 1/rho; % Diagonal and with the same value: scalar
        phi = zeros(length(x),length(t)); % No initial damage
        % Initial condition
        u_n = zeros(2*length(x),length(t)); % 2*N x n  2D matrix
        v_n = zeros(2*length(x),2); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.ffz
        if size(body_force,2) == 1
            % Constant body force
            bn = [body_force body_force]; % n n+1
            flag_bf_constant = true;
        else
            if size(body_force,2) ~= length(t)
                bfnew = zeros(2*length(x),length(t));
                for iii = 1:size(body_force,1)
                    bfnew(iii,:) = interp1(0:size(body_force,2)-1,body_force(iii,:),(0:length(t)-1)*(size(body_force,2)-1)/(length(t)-1));
                end
                body_force = bfnew;
            end
            flag_bf_constant = false;
        end
        n_sample = [1 data_dump:data_dump:length(t)];
        if n_sample(end)~=length(t)
            n_sample = [n_sample length(t)];
        end
        t_s = t(n_sample);
        energy.W = zeros(length(x),length(n_sample)); % No initial  deformation hence no initial strain energy density
        energy.KE = zeros(length(x),length(n_sample));
        energy.EW = zeros(length(x),length(n_sample));
        energy_ext_var = zeros(length(x),1); % External energy that takes into account variable body force (velocity constraints)
        fn = zeros(2*length(x),1); % Initial force
        u_const = zeros(length(v_n)-(ndof),1); % Constraint nodes
        % Temporary variables
        history_S = model.history.S; %
         if model.dilatation 
             history_T = model.history.theta;
         end
        fn_temp = zeros(size(x));
        phi_temp = zeros(length(x),1);
        
        % {Recoverying temporary files}
        n_initial = 1;
        if exist('tempsim.mat','file')
            load('tempsim.mat');
            n_initial = n;
            disp('Found a partial simulation. Continuing it...')
        end
        timerVal = tic; % Initiating another tic
        for n = n_initial:length(t)-1
            % Instatiate body force
            if ~flag_bf_constant
                bn = body_force(:,n:n+1); % Increment b
            end
            %% ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % #### Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof,1)); % V(n+1/2)
            % #### Step 2 - Update displacement
            u_n(1:ndof,n+1) = u_n(1:ndof,n) + dt*v_n(1:ndof,2); % u(n+1) - %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
             % ----- Solving for the displacement constraint nodes ----
            if ~isempty(u_const)
                if sum(bc_set(:,3) ~= 0 & bc_set(:,2) ~= 0)
                    error('The boundary condition matrix bc_set has prescribed both displacement and velocity for the same node.')
                end
                u_const(bc_set(:,3) == 0) = bc_set(bc_set(:,3) == 0,2); % Defining the displacements for the nodes with no velocity
                u_const = u_const + bc_set(:,3)*dt; % Updating the velocity constraint nodes
                u_n(ndof+1:end,n+1) = u_const;
            end
            u1 = u_n(:,n);
            u2 = u_n(:,n+1); % Vector of displacement that will be used to come back to 2D
            
            % ---- {Evaluating dilatation} ----
            theta = zeros(length(x),1); % Preallocate theta
            damage.phi = phi(:,n); % Accessing current damage situation
            if model.dilatation
                if b_parll
                    parfor ii = 1:length(x)
                        [theta(ii),history_T(ii)] = model.dilatationEval(x,u2,familyMat(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,damage,history_T(ii));
                    end
                else
                    [theta,history_T] = model.dilatationEval(x,u2,familyMat,partialAreas,surfaceCorrection,[],idb,par_omega,damage,history_S,history_T);
                end
                %model.history.theta = history_tempT; % Assigning up-to-date history variable
            end
            
            % ####### Step 3 - Update velocity 
            % {Evaluate f[n+1]}
            fn = zeros(2*length(x),1); % Instatiate force vector
            energy_pot = zeros(length(x),1); % Pre-allocate potential energy
            energy_ext = zeros(length(x),1); % Pre-allocate external force 
            b_Weval = rem(n+1,data_dump) == 0 || n+1 == length(t); % Deciding when to evaluate the energy
            if b_parll
                parfor ii = 1:length(x)
                   [fn_temp(ii,:),history_S(ii,:),phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,phi(:,n),history_S(ii,:),history_T,dt,A,body_force,theta,b_Weval,bc_set);
                end
            else 
                for ii = 1:length(x)
                   [fn_temp(ii,:),history_S(ii,:),phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,phi(:,n),history_S(ii,:),history_T,dt,A,body_force,theta,b_Weval,bc_set);
                end
            end
            % Converting the temporary variables
            for ii = 1:length(x)
                fn(dof_vec(ii,:)) = fn_temp(ii,:)';
            end
            phi(:,n+1) = phi_temp;
            % Evaluate V(n+1)
            v_n(1:ndof,1) = v_n(1:ndof,2) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof,1)); % V(n+1) is stored in the next V(n) using f n+1 and b n+1
            
            %% Evaluating energy
            if b_Weval
                index_s = n_sample == n+1;
                % {Potential energy}
                energy.W(:,index_s) = energy_pot;
                % {External work}
                BBN = bn(:,1);
                if flag_bf_constant
                    % Constant body force
                    energy_ext = dot(u2(dof_vec)',BBN(dof_vec)')'.*A;
                else
                    du = u2(dof_vec) - u1(dof_vec);
                    energy_ext_var = energy_ext_var+dot(du',BBN(dof_vec)')'.*A;
                end
                % {External work realized by the velocity constraint}
                if ~isempty(bc_set)
                    vel_dof = idb(bc_set(bc_set(:,3)~=0,1));
                    v = bc_set(bc_set(:,3)~=0,3);
                else
                    vel_dof = [];
                    v = [];
                end
                bf = -fn(vel_dof);
                du = v*dt;
                add_ext = bf.*du;
                for gg = 1:length(vel_dof)
                    ind_vel = find(dof_vec == vel_dof(gg));
                    ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
                    energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
                end
                %energy_ext_var(vel_dof) = energy_ext_var(vel_dof) + bf.*du*V;
                energy.EW(:,index_s) = energy_ext + energy_ext_var;%+energy.EW(:,n);
                % Kinectic energy
                for kk = 1:length(x)
                    dofk = dof_vec(kk,:);
                    energy.KE(kk,index_s) =  1/2*rho*norm(v_n(dofk,1))^2.*A(kk);
                end
            else
                 % {External incremental work only}
                BBN = bn(:,1);
                if~ flag_bf_constant
                    du = u2(dof_vec) - u1(dof_vec);
                    energy_ext_var = energy_ext_var+dot(du',BBN(dof_vec)')'.*A;
                end
                % {External work realized by the velocity constraint}
                if ~isempty(bc_set)
                    vel_dof = idb(bc_set(bc_set(:,3)~=0,1));
                    v = bc_set(bc_set(:,3)~=0,3);
                else
                    vel_dof = [];
                    v = [];
                end
                bf = -fn(vel_dof);
                du = v*dt;
                add_ext = bf.*du;
                for gg = 1:length(vel_dof)
                    ind_vel = find(dof_vec == vel_dof(gg));
                    ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
                    energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
                end
            end
            
            %% ############ COUNTING THE PROCESSING TIME #############
            time_up = toc(timerVal);
            clc
            disp("Time = " + num2str(t(n)) + " secs. Percentage of the process: " + num2str(n/(length(t)-1)*100) + "%. ETA: "+ num2str(time_up/n*(length(t)-n)))         
        end
        % Sampling the results
        [u_n] = sampling(u_n,t,t_s);
        [phi] = sampling(phi,t,t_s);
end
%%
function dt_crit = criticalTimeStep(x,family,partialAreas,par_omega,c,rho,model)
    dt = zeros(length(x),1);
    for ii = 1:length(x)
        if model.name == "PMB" || model.name == "PMB DTT" || model.name == "LBB"
            familyOfI = family(ii,family(ii,:)~=0);
            den = 0;
            for jj = familyOfI
                neigh_index = family(ii,:) == jj;
                xi = x(jj,:) - x(ii,:);
                norma = norm(xi);
                C = c(1)*influenceFunction(norma,par_omega)*norma/norma^3*...
                [xi(1)^2 xi(1)*xi(2); xi(1)*xi(2) xi(2)^2]; % PMB model only
                den = den + norm(C)*partialAreas(neigh_index);
            end
            dt(ii) = sqrt(2*rho/den);
        else
            break;
        end
    end
    dt_crit = min(dt); % Critical time step
end
%%
function [f_i,history_S_up,phi_up,energy_pot] = parFor_loop(x,u_n,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,model,damage,phi,history_S,history_T,dt,A,body_force,theta,b_Weval,bc_set)
   % Loop on the nodes
   family = familyMat(ii,familyMat(ii,:)~=0);
   neig_index = 1:length(family);
   jj = family(neig_index);
   % Loop on their neighbourhood
   history_S_up = history_S; 
   noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
   if model.dilatation
      [fij,history_S_up(neig_index),mu_j] = model.T(x,u_n,theta,ii,jj,dof_vec,par_omega,[ ],damage,history_S(neig_index),history_T,noFail);
   else
      [fij,history_S_up(neig_index),mu_j] = model.T(x,u_n,ii,jj,dof_vec,par_omega,[ ],damage,history_S(neig_index),dt,noFail);
   end
   Vj = partialAreas(ii,neig_index)';
   lambda = surfaceCorrection(ii,neig_index)';
   f_i = sum(fij.*Vj.*lambda);
   % Damage index
   areaTot = sum(Vj);
   partialDamage = sum(mu_j.*Vj);
   phi_up = 1 - partialDamage/areaTot;
   if b_Weval
       if ~model.dilatation
           % Strain energy
           W = model.strainEnergyDensity(x,u_n,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,damage,history_S_up(neig_index));
       else
           W = model.strainEnergyDensity(x,u_n,theta,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,damage,history_S_up(neig_index),history_T); % neig_index == length(family)
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