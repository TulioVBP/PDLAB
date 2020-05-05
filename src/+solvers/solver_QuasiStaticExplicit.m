function [t_s,u_load,u_n,n,phi,energy,history,time_up] = solver_QuasiStaticExplicit(x,t,idb,body_force,...
   bc_set,familyMat,A,partialAreas,surfaceCorrection,T,c,model,par_omega,history,noFailZone,damage,b_parll,beta,n_load,data_dump)
% Explicit time solver with dynamic relaxation for quasi-static experiments
   if nargin < 18 % No data dump
        data_dump = 1;
    elseif nargin < 20
        data_dump = 1;
   end
   
   if isempty(t)
       t = 1:10000;
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
    damage.checkCrack = zeros(size(history.S));
    for ii = 1:size(history.S)
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
    phi = zeros(length(x),length(t)); % No initial damage
    % Initial condition
    u_n = zeros(2*length(x),length(t)); % 2*N x n  2D matrix
    v_n = zeros(2*length(x),2); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.
    
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
    
     
    ndof = 2*length(x) - size(bc_set,1); % Number of degree of freedoms
    
    % Defining n_sample and initial energy variables
    n_sample = [1 data_dump:data_dump:length(t)];
    if n_sample(end)~=length(t)
        n_sample = [n_sample length(t)];
    end
    t_s = t(n_sample); % Sample time
    energy.W = zeros(length(x),length(n_sample)); % No initial  deformation hence no initial strain energy density
    energy.KE = zeros(length(x),length(n_sample));
    energy.EW = zeros(length(x),length(n_sample));
    energy_ext_var = zeros(length(x),1); % External energy that takes into account variable body force (velocity constraints)
    fn = zeros(2*length(x),1); % Initial force
    u_const = zeros(length(v_n)-(ndof),1); % Constraint nodes
    
    % Temporary variables
    history_tempS = history.S;
    if model.dilatation 
        history_tempT = history.theta;
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
    u_load = zeros(2*length(x),n_load);
    for kk = 1:n_load
        %% #########  LOAD LOOP ##########
        bc_set_part = bc_set;
        bc_set_part(:,2) = bc_set(:,2)*kk/n_load; % Upload boundary condition
        body_force_part(:,1) = body_force*kk/n_load; % bn
        crit_var = zeros(1,length(t)-1);
        % Density vector
        dt = abs(t(2)-t(1));
        lambda = evaluateLambda(x,[],ndof,idb,familyMat,partialAreas,surfaceCorrection,ones(2*length(x)),par_omega,c,model,damage,history,dt);
        lambda_inv = 1./lambda;
        C = 0;
        F = zeros(2*length(x),2);
        for n = n_initial:length(t)-1
            % Instatiate body force
            bn = body_force_part; % Increment b
            %% ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % #### Step 0 - Acceleration at a_n
            an =  lambda_inv .* (fn(1:ndof) + bn(1:ndof) - C * lambda.* v_n(1:ndof,1));
            % #### Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt/2 * an; % V(n+1/2)
            % #### Step 2 - Update displacement
            u_n(1:ndof,n+1) = u_n(1:ndof,n) + dt*v_n(1:ndof,2); % u(n+1) - %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
             % ----- Solving for the displacement constraint nodes ----
            if ~isempty(u_const)
                u_const = bc_set_part(:,2); % Defining the displacements for the nodes with no velocity
                %u_const = u_const + bc_set(:,3)*dt; % Updating the velocity constraint nodes
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
                        [theta(ii),history_tempT(ii)] = dilatation(x,u2,familyMat(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,c,model,damage,history,dt);
                    end
                else
                    [theta,history_tempT] = dilatation(x,u2,familyMat,partialAreas,surfaceCorrection,[],idb,par_omega,c,model,damage,history,dt);
                end
                history.theta = history_tempT; % Assigning up-to-date history variable
            end

            % ####### Step 3 - Update velocity 
            % Storing fn
            F(:,1) = fn + bn; % F(n)
            % {Evaluate f[n+1]}
            fn = zeros(2*length(x),1); % Instatiate force vector
            energy_pot = zeros(length(x),1); % Pre-allocate potential energy
            energy_ext = zeros(length(x),1); % Pre-allocate external force 
            b_Weval = rem(n+1,data_dump) == 0 || n+1 == length(t); % Deciding when to evaluate the energy
            if b_parll
                parfor ii = 1:length(x)
                   [fn_temp(ii,:),history_tempS(ii,:),phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi(:,n),dt,history,T,A,body_force,theta,b_Weval,bc_set);
                end
            else 
                for ii = 1:length(x)
                   [fn_temp(ii,:),history_tempS(ii,:),phi_temp(ii),energy_pot(ii)] = parFor_loop(x,u2,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi(:,n),dt,history,T,A,body_force,theta,b_Weval,bc_set);
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
            C = evaluateDamping(lambda,u_n(:,n+1),ndof,dt,F,v_n(:,2));
            
%             %% Evaluating energy
%             if b_Weval
%                 index_s = n_sample == n+1;
%                 % {Potential energy}
%                 energy.W(:,index_s) = energy_pot;
%                 % {External work}
%                 if flag_bf_constant
%                     % Constant body force
%                     energy_ext = dot(u2(dof_vec)',bn(dof_vec)')'.*A;
%                 else
%                     du = u2(dof_vec) - u1(dof_vec);
%                     energy_ext_var = energy_ext_var+dot(du',bn(dof_vec)')'.*A;
%                 end
%                 % {External work realized by the velocity constraint}
%                 if ~isempty(bc_set)
%                     vel_dof = idb(bc_set(bc_set(:,3)~=0,1));
%                     v = bc_set(bc_set(:,3)~=0,3);
%                 else
%                     vel_dof = [];
%                     v = [];
%                 end
%                 bf = -fn(vel_dof);
%                 du = v*dt;
%                 add_ext = bf.*du;
%                 for gg = 1:length(vel_dof)
%                     ind_vel = find(dof_vec == vel_dof(gg));
%                     ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
%                     energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
%                 end
%                 %energy_ext_var(vel_dof) = energy_ext_var(vel_dof) + bf.*du*V;
%                 energy.EW(:,index_s) = energy_ext + energy_ext_var;%+energy.EW(:,n);
%                 % Kinectic energy
%                 for kk = 1:length(x)
%                     dofk = dof_vec(kk,:);
%                     energy.KE(kk,index_s) =  1/2*rho*norm(v_n(dofk,1))^2.*A(kk);
%                 end
%             else
%                  % {External incremental work only}
%                 if~ flag_bf_constant
%                     du = u2(dof_vec) - u1(dof_vec);
%                     energy_ext_var = energy_ext_var+dot(du',bn(dof_vec)')'.*A;
%                 end
%                 % {External work realized by the velocity constraint}
%                 if ~isempty(bc_set)
%                     vel_dof = idb(bc_set(bc_set(:,3)~=0,1));
%                     v = bc_set(bc_set(:,3)~=0,3);
%                 else
%                     vel_dof = [];
%                     v = [];
%                 end
%                 bf = -fn(vel_dof);
%                 du = v*dt;
%                 add_ext = bf.*du;
%                 for gg = 1:length(vel_dof)
%                     ind_vel = find(dof_vec == vel_dof(gg));
%                     ind_vel  = ind_vel - (ind_vel > size(dof_vec,1))*size(dof_vec,1);
%                     energy_ext_var(ind_vel) = energy_ext_var(ind_vel)+ add_ext(gg).*A(ind_vel);
%                 end
%             end

            %% ############ Verifying convergence ##########
            PHI = fn(1:ndof) + bn(1:ndof);
            crit_var(n) = norm(PHI)/norm(bn(1:ndof));
            if crit_var(n) < beta
                disp("Convergence achieved for the load step " + num2str(kk) + " ...")
                break
            else
                 %% ############ COUNTING THE PROCESSING TIME #############
                time_up = toc(timerVal);
                if n >1
                   clc
                end
                %lineLength = fprintf("Load step " + num2str(kk) + " out of " + num2str(n_load) +": Time step "+ num2str(n) + "%. ETA: "+ num2str(time_up/n*(length(t)-n)));
                disp("Load step " + num2str(kk) + " out of " + num2str(n_load) +": Time step "+ num2str(n) + ". Criterion: " + crit_var(n) + ". ETA: "+ num2str(time_up/n*(length(t)-n)))        
            end
            if n == length(t)-1
                figure
                plot(1:length(t)-1,crit_var)
                hold on
                plot(1:length(t),beta*ones(length(t),1))
                grid on
                xlabel('Time step')
                ylabel('\psi')
                pause
                
                %error('Failure to converge with given time steps.')
            end
            %% Verifying rate of convergence
            if n > 2
                if (abs(crit_var(n) - crit_var(n-1)))/crit_var(n-1) < 10^-10
                    error("Convergence rate is too slow and convergence criteria was not met. Conv. crit. = " + num2str(crit_var(n)));
                end
            end
        end
        u_load(:,kk) = u_n(:,n);
    end
end

%%
function [f_i,history_upS,phi_up,energy_pot] = parFor_loop(x,u_n,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi,dt,history,T,A,body_force,theta,b_Weval,bc_set)
   %dofi = dof_vec(ii,:);
   % Loop on the nodes
   %areaTot = 0; partialDamage = 0; % Instatiate for damage index
   family = familyMat(ii,familyMat(ii,:)~=0);
   %f_i = 0;
   history_upS = history.S(ii,:);
   %damage.phi = phi(ii);
   neig_index = 1:length(family);
   jj = family(neig_index);
   % Loop on their neighbourhood
   noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
   if model.dilatation
      %[fij,history_upS(neig_index),mu_j] = T(x,u_n,ii,dof_vec,familyMat,partialAreas,neig_index,par_omega,c,model,[ ],damage,dt,history.S(ii,neig_index),history.theta,noFail);
      [fij,history_upS(neig_index),mu_j] = T(x,u_n,theta,ii,jj,dof_vec,par_omega,c,model,[ ],damage,dt,history.S(ii,neig_index),history.theta,noFail);
   else
      [fij,history_upS(neig_index),mu_j] = T(x,u_n,ii,jj,dof_vec,par_omega,c,model,[ ],damage,dt,history.S(ii,neig_index),noFail);
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
           W = strainEnergyDensity(x,u_n,[],familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,c,model,damage,history_upS(1,neig_index),[]);
       else
           W = strainEnergyDensity(x,u_n,theta,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,c,model,damage,history_upS(neig_index),history.theta); % neig_index == length(family)
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

function lambda = evaluateLambda(x,u,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history,dt)
    % 1 - Define stiffness matrix
    K = analyticalStiffnessMatrix(x,u,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history);
    % 2 - Define lambda
    lambda = diag(eye(ndof));
    for ii = 1:ndof
        lambda(ii) = lambda(ii) * (1/4*dt^2*sum(abs(K(ii,1:ndof))));
    end
end

function C = evaluateDamping(lambda,u,ndof,dt,F,v)
    % 4 - Define K diagonal
    vv = v(1:ndof);
    Kn = -(F(1:ndof,2) - F(1:ndof,1))./lambda./(dt*vv);
    Kn (vv < 1e-14 & vv > -1e-14) = 0;
    Kn = diag(Kn);
    % 5 - Evaluate the damping
    if u(1:ndof)'*Kn*u(1:ndof) < 0
        disp("Negative Kn");
    end
    C = real(2*sqrt(u(1:ndof)'*Kn*u(1:ndof)/(u(1:ndof)'*u(1:ndof))));
    disp("C is "+num2str(C));
    if C < 0 
        C = 0;
    end
end