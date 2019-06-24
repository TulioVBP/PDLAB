% Explicit time - dynamic solver

function [u_n,phi,energy] = solver_DynamicExplicit(x,t,idb,body_force,bc_set,familyMat,partialAreas,surfaceCorrection,T,c,rho,model,par_omega,history,noFailZone,damage,b_parll)
    ndof = 2*length(x) - length(bc_set);
    %% Create volume
    h = norm(x(1,:)-x(2,:));
    V = h^2; % homogeneous volume
    %% Verify if dt is small enough
    dt = abs(t(2) - t(1));
    dt_crit = criticalTimeStep(x,familyMat,partialAreas,par_omega,c,rho,model);
    dt_ratio = dt/dt_crit;
    if dt_ratio < 1
        disp("Given time-step " + num2str(dt) +" sec is less than the critical time-step, " + num2str(dt_crit) + ...
            " sec, and the safety factor is " + num2str(dt_ratio) + ". The simulation should converge.")
    else 
       disp("Given time-step " + num2str(dt) + " sec is greater than the critical time-step, " + num2str(dt_crit) + ...
            " sec, and the safety factor is " + num2str(dt_ratio)+ ". The simulation shall explode.")
    end
    %% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
    %% INITIALIZE SIMULATION MATRICES
        if b_parll
            poolobj = parpool(2);
        end
        Minv = 1/rho; % Diagonal and with the same value: scalar
        phi = zeros(length(x),length(t)); % No initial damage
        W = zeros(length(x),length(t)); % No initial deformation
        % Initial condition
        u_n = zeros(2*length(x),length(t)); % 2*N x n  2D matrix
        v_n = zeros(2*length(x),2); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.ffz
        bn = body_force;%bodyForce(x,stresses,m,h,A);
        energy.W = zeros(length(x),length(t));
        energy.KE = zeros(length(x),length(t));
        energy.EW = zeros(length(x),length(t));
        fn = zeros(2*length(x),1); % Initial force
        u_const = zeros(length(v_n)-(ndof),1); % Constraint nodes
        % Temporary variables
        history_temp = history;
        fn_temp = zeros(size(x));
        phi_temp = zeros(length(x),1);
        % {Recoverying temporary files}
        n_initial = 1;
        if exist('tempsim.mat','file')
            load('tempsim.mat');
            n_initial = n;
            disp('Found a partial simulation. Continuing it...')
        end
        for n = n_initial:length(t)-1
            %% ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % #### Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof)); % V(n+1/2)
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
            % ### Step 3 - Update velocity
            % Evaluate f[n+1]
            fn = zeros(2*length(x),1); % Instatiate force vector
            if b_parll
                parfor ii = 1:length(x)
                   [fn_temp(ii,:),history_temp(ii,:,:),phi_temp(ii),energy_pot(ii),energy_ext(ii)] = parFor_loop(x,u_n(:,n+1),dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi(:,n),dt,noFailZone,history,T,W,V,body_force);
                end
            else 
                for ii = 1:length(x)
                   [fn_temp(ii,:),history_temp(ii,:,:),phi_temp(ii),energy_pot(ii),energy_ext(ii)] = parFor_loop(x,u_n(:,n+1),dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi(:,n),dt,noFailZone,history,T,W,V,body_force);
                end
            end
            % Converting the temporary variables
            for ii = 1:length(x)
                fn(dof_vec(ii,:)) = fn_temp(ii,:)';
            end
            history = history_temp;
            phi(:,n+1) = phi_temp;
            energy.W(:,n+1) = energy_pot;
            energy.EW(:,n+1) = energy_ext;
            % Evaluate V(n+1)
            v_n(1:ndof,1) = v_n(1:ndof,2) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof)); % V(n+1) is stored in the next V(n)
            % Kinectic energy
            for kk = 1:length(x)
                dofk = dof_vec(kk,:);
                energy.KE(kk,n+1) =  1/2*rho*norm(v_n(dofk,1))^2*V;
            end
            
            % ############ COUNTING THE PROCESSING TIME #############
            disp("Time = " + num2str(t(n)) + " secs. Percentage of the process: " + num2str(n/(length(t)-1)*100) + "%")
            %pause
            % {Checking for running out time}
            horario = clock;
            if horario(4) >= 23 && horario(5) >40
                filename = strcat('tempsim.mat');
                save(filename,'x','idb','u_n','phi','energy','n');
                if b_parll
                    delete(poolobj);
                end
                break;
            end
                
        end
        if b_parll
            delete(poolobj);
        end
end
%%
function dt_crit = criticalTimeStep(x,family,partialAreas,par_omega,c,rho,model)
    dt = zeros(length(x),1);
    for ii = 1:length(x)
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
    end
    dt_crit = min(dt); % Critical time step
end
%%
function [f_i,history_up,phi_up,energy_pot,energy_ext] = parFor_loop(x,u_n,dof_vec,idb,ii,familyMat,partialAreas,surfaceCorrection,par_omega,c,model,damage,phi,dt,noFailZone,history,T,W,V,body_force)
   dofi = dof_vec(ii,:);
   % Loop on the nodes
   areaTot = 0; partialDamage = 0; % Instatiate for damage index
   family = familyMat(ii,familyMat(ii,:)~=0);
   f_i = 0;
   history_up = history(ii,:,:);
   damage.phi = phi(ii);
   for neig_index = 1:length(family)
       jj = family(neig_index);
       % Loop on their neighbourhood
       noFail = noFailZone(ii) || noFailZone(jj); % True if node ii or jj is in the no fail zone
       if model.dilatation
          [fij,history_up(1,neig_index,:),mu_j] = T(x,u_n,ii,dof_vec,familyMat,partialAreas,neig_index,par_omega,c,model,[ ],damage,dt,history(ii,neig_index,:),noFail);
       else
          [fij,history_up(1,neig_index,:),mu_j] = T(x,u_n,ii,jj,dof_vec,par_omega,c,model,[ ],damage,dt,history(ii,neig_index),noFail);
       end
       Vj = partialAreas(ii,neig_index);
       lambda = surfaceCorrection(ii,neig_index);
       f_i = f_i + (fij)*Vj*lambda;
       % Damage index
       areaTot = areaTot + Vj;
       partialDamage = partialDamage + mu_j*Vj;
        %phi(ii,n+1) = phi(ii,n+1) - (damageIndex(x,u_n(:,n+1),familyMat(ii,neig_index),partialAreas(ii,neig_index),ii,idb,noFailZone)-1); % Damage index
       if ~model.dilatation
           % Strain energy
           W(ii) = W(ii) + strainEnergyDensity(x,u_n,familyMat(ii,neig_index),partialAreas(ii,neig_index),surfaceCorrection(ii,neig_index),ii,idb,par_omega,c,model,damage,history_up(1,neig_index));
           % 1/2 factor not truly understood.
       end
   end
   if model.dilatation
       % Strain energy
       W(ii) = strainEnergyDensity(x,u_n,familyMat(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,c,model,damage,history_up(1,:,:));
   end
   phi_up = 1 - partialDamage/areaTot;
   % External work
   energy_ext = dot(u_n(dofi),body_force(dofi))*V;
   % Stored strain energy
   energy_pot = W(ii)*V;
   % Kinectic energy - 
   %energy.KE(ii,n+1) =  1/2*rho*norm(v_n(dofi,2))^2*V;
end
