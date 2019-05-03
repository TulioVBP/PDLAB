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
        Minv = 1/rho/V; % Diagonal and with the same value: scalar
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
        for n = 1:length(t)-1
            %% TIME LOOP
            fn = zeros(2*length(x),1);
            for ii = 1:length(x)
                dofi = dof_vec(ii,:);
                % Loop on the nodes
               family = familyMat(ii,familyMat(ii,:)~=0);
               for jj = family
                   % Loop on their neighbourhood
                   noFail = noFailZone(ii) || noFailZone(jj); % True if node ii or jj is in the no fail zone
                   neig_index = find(familyMat(ii,:) == jj);
                   %[fij,history] = T(x(ii,:),x(jj,:),u_n(dofi,n)',u_n(dofj,n)',dt,history(ii,neig_index,:),noFail)
                   [fij,history(ii,neig_index)] = T(x,u_n(:,n),ii,dof_vec,familyMat,neig_index,dt,history(ii,neig_index),noFail);
                   Vj = partialAreas(ii,neig_index);
                   fn(dofi) = fn(dofi) + (fij')*Vj;
                   % Damage index
                    phi(ii,n) = phi(ii,n) + damageIndex(x,u_n(:,n),familyMat(ii,neig_index),partialAreas(ii,neig_index),ii,idb,noFailZone); % Damage index
                   % Strain energy
                    W(ii,n) = W(ii,n) + strainEnergyDensity(x,u_n(:,n),familyMat(ii,neig_index),partialAreas(ii,neig_index),ii,history(ii,neig_index),idb);
               end
               % Kinectic energy
               energy.KE(ii,n) =  1/2*rho*norm(v_n(dofi,1))^2*V;
               % External work
               energy.EW(ii,n) = dot(u_n(dofi,n),b(dofi))*V;
               % Stored strain energy
               energy.W(ii,n) = W(ii,n)*V;
               %phi(ii,n) = damageIndex(x,u_n(:,(2*n-1):(2*n)),family(ii,:),partialAreas(ii,:),ii,crackIn); % Damage index
            end
            % ############ VELOCITY VERLET ALGORITHM ###############
            % ---- Solving for the dof ----
            % Step 1 - Midway velocity
            v_n(1:ndof,2) = v_n(1:ndof,1) + dt*Minv*(fn(1:ndof)*V + bn(1:ndof)*V); % V(n+1/2)
            %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
            u_n(1:ndof,n+1) = u_n(1:ndof,n) + dt*v_n(1:ndof,1); % u(n+1)
            v_n(1:ndof,1) = v_n(1:ndof,2); %+ dt*Minv*(fn(1:ndof)*V + bn(1:ndof)*V); % V(n+1) is stored in the next V(n)
            % ----- Solving for the constraint nodes ----
            u_const = zeros(length(v_n)-(ndof),1);
            if ~isempty(u_const)
                u_const(bc_set(:,3) == 0) = bc_set(bc_set(:,3) == 0,2); % Defining the displacements for the nodes with no velocity
                u_const = u_const + bc_set(:,3)*dt;
                u_n(ndof+1:end,n+1) = u_n(ndof+1:end,n) + u_const;
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

