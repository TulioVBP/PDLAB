% Explicit time - dynamic solver

function [u_n,phi,energy] = solver_DynamicExplicit(x,t,idb,b,bc_set,familyMat,partialAreas,T,history,noFailZone)
    global rho
    ndof = 2*length(x) - length(bc_set);
    %% Create volume
    h = norm(x(1,:)-x(2,:));
    V = h^2; % homogeneous volume
    %% Verify if dt is small enough
    dt = abs(t(2) - t(1));
    dt_crit = criticalTimeStep(x,familyMat,partialAreas);
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
        Minv = 1/rho; % Diagonal and with the same value: scalar
        %S_max = zeros(length(x),history);
        phi = zeros(length(x),length(t)); % No initial damage
        W = zeros(length(x),length(t)); % No initial deformation
        % Initial condition
        u_n = zeros(2*length(x),length(t)); % 2*N x n  2D matrix
        v_n = zeros(2*length(x),2); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.ffz
        bn = b;%bodyForce(x,stresses,m,h,A);
        energy.W = zeros(length(x),length(t));
        energy.KE = zeros(length(x),length(t));
        energy.EW = zeros(length(x),length(t));
        fn = zeros(2*length(x),1); % Initial force
        for n = 1:length(t)-1
            %% TIME LOOP
            % ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % #### Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof)); % V(n+1/2)
            % #### Step 2 - Update displacement
            u_n(1:ndof,n+1) = u_n(1:ndof,n) + dt*v_n(1:ndof,2); % u(n+1) - %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
             % ----- Solving for the constraint nodes ----
            u_const = zeros(length(v_n)-(ndof),1);
            if ~isempty(u_const)
                u_const(bc_set(:,3) == 0) = bc_set(bc_set(:,3) == 0,2); % Defining the displacements for the nodes with no velocity
                u_const = u_const + bc_set(:,3)*dt;
                u_n(ndof+1:end,n+1) = u_n(ndof+1:end,n) + u_const;
            end
            % ### Step 3 - Update velocity
            % Evaluate f[n+1]
            fn = zeros(2*length(x),1); % Instatiate force vector
            for ii = 1:length(x)
                dofi = dof_vec(ii,:);
                % Loop on the nodes
               family = familyMat(ii,familyMat(ii,:)~=0);
               for jj = family
                   % Loop on their neighbourhood
                   noFail = noFailZone(ii) || noFailZone(jj); % True if node ii or jj is in the no fail zone
                   neig_index = find(familyMat(ii,:) == jj);
                   %[fij,history] = T(x(ii,:),x(jj,:),u_n(dofi,n)',u_n(dofj,n)',dt,history(ii,neig_index,:),noFail)
                   [fij,history(ii,neig_index)] = T(x,u_n(:,n+1),ii,dof_vec,familyMat,partialAreas,neig_index,dt,history(ii,neig_index),noFail);
                   Vj = partialAreas(ii,neig_index);
                   fn(dofi) = fn(dofi) + (fij')*Vj;
                   % Damage index
                    phi(ii,n+1) = phi(ii,n+1) + damageIndex(x,u_n(:,n+1),familyMat(ii,neig_index),partialAreas(ii,neig_index),ii,idb,noFailZone); % Damage index
                   % Strain energy
                    W(ii,n+1) = W(ii,n+1) + strainEnergyDensity(x,u_n(:,n+1),familyMat(ii,neig_index),partialAreas(ii,neig_index),ii,history(ii,neig_index),idb);
               end
               % External work
               energy.EW(ii,n+1) = dot(u_n(dofi,n+1),b(dofi))*V;
               % Stored strain energy
               energy.W(ii,n+1) = W(ii,n+1)*V;
               % Kinectic energy - 
               %energy.KE(ii,n+1) =  1/2*rho*norm(v_n(dofi,2))^2*V;
            end
            v_n(1:ndof,1) = v_n(1:ndof,2) + dt/2*Minv*(fn(1:ndof) + bn(1:ndof)); % V(n+1) is stored in the next V(n)
            % Kinectic energy
            for kk = 1:length(x)
                dofk = dof_vec(kk,:);
                energy.KE(:,n+1) =  1/2*rho*norm(v_n(dof_vec,1))^2*V;
            end
            % ############ COUNTING THE PROCESSING TIME #############
            disp("Time = " + num2str(t(n)) + " secs. Percentage of the process: " + num2str(n/(length(t)-1)*100) + "%")
            %pause
        end
end

function dt_crit = criticalTimeStep(x,family,partialAreas)
    global rho c1 omega horizon
    dt = zeros(length(x),1);
    for ii = 1:length(x)
        familyOfI = family(ii,family(ii,:)~=0);
        den = 0;
        for jj = familyOfI
            neigh_index = family(ii,:) == jj;
            xi = x(jj,:) - x(ii,:);
            norma = norm(xi);
            C = c1*influenceFunction(norma,horizon,omega)*norma/norma^3*...
            [xi(1)^2 xi(1)*xi(2); xi(1)*xi(2) xi(2)^2]; % PMB model only
            den = den + norm(C)*partialAreas(neigh_index);
        end
        dt(ii) = sqrt(2*rho/den);
    end
    dt_crit = min(dt); % Critical time step
end

