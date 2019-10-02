%QUASI-STATIC SOLVER

function [un,r,energy,phi] = solver_QuasiStatic(x,n_tot,idb,b,bc_set,family,partialAreas,surfaceCorrection,T,c,model,par_omega,ndof,V,damage,history,noFailZone)
N = length(idb);
%% Damage variables 
%{Defining cracking trespassing matrix}
    %if damage.damageOn
        crackSegments = size(damage.crackIn,1); % At least 2
        damage.checkCrack = zeros(size(history.S));
        for ii = 1:size(history.S)
            x_j = x(family(ii,family(ii,:)~=0),:);
            check = zeros(size(x_j,1),crackSegments-1);
            for kk = 1:crackSegments-1
                for jj = 1:size(x_j,1)
                [~,check(jj,kk)] = checkBondCrack(x(ii,:),x_j(jj,:),damage.crackIn(kk,:),damage.crackIn(kk+1,:));
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
    %end
    % {No fail to damage variable}
    damage.noFail = noFailZone;   
    phi = zeros(length(x),n_tot); % Damage index
    

%% Step 1 - Initialization
un = zeros(N,n_tot); % N= 2*nn
energy.W = zeros(length(x),n_tot);
energy.KE = zeros(length(x),n_tot);
energy.EW = zeros(length(x),n_tot);
% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
for n = 1:n_tot
    %% Step 2 - Update the load step n <- n + 1 and pseudo-time t. Update the boundary conditions.
    bn = b*(n/n_tot); % Partial load
    %% Step 3 - Evaluate the residual vector, r, and residual r. Determine the convergence criterion
    %          for the load step.
    epsilon = 10^-4;
    damage.phi = phi(:,n); % Accessing current damage situation
    [r_vec,history,phi(:,n)] = getForce(x,un(:,n),T,bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V,par_omega,c,model,damage,history); % Update to include arbitrary displacement kinematic conditions
    r_max = epsilon*max(norm(bn*V,Inf),norm(r_vec-bn*V,Inf)); % Normalizing the maximum residual
    %% Step 4 - Assign an initial guess to the trial displacement utrial (for example, utrial = un).
    u_trial = un(:,n);
    %% Step 5 - Apply Newton's method to minimize the residual.
    r = norm(r_vec,Inf);
    % -------------------- Newton's method ----------------
    if ~model.linearity
        % Suitable for non-linear models
        while r > r_max
            % {Damage}
            damage.phi = phi(:,n); % Accessing current damage situation
            if ~model.stiffnessAnal 
                K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,surfaceCorrection,T,ndof,par_omega,c,model,damage,history);
            else
                K = analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history);
            end
            disp('Stiffness matrix done.')
            du = -K\r_vec;
            disp('Incremental solution found.')
            u_trial = u_trial + du;
            [r_vec,history,phi(:,n)] = getForce(x,u_trial,T,bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V,par_omega,c,model,damage,history); % Update to include arbitrary displacement kinematic conditions
            r = norm(r_vec,Inf);
            %r_max = epsilon*max(norm(bn*V,Inf),norm(r_vec-bn*V,Inf));
            disp("Current residual equal to "+ num2str(r) + ". Maximum residual to be " + num2str(r_max))
        end
        disp("Solution found for the step " + int2str(n) + " out of " + int2str(n_tot))
        un(:,n) = u_trial;
        if n < n_tot
            un(:,n+1) = u_trial; % For the next iteration
        end
    else
        % Linear model
        % {Damage}
        damage.phi = phi(:,n); % Accessing current damage situation
        if n < 2 % If the model is linear, there is no need to find the matrix more than once
            if ~model.stiffnessAnal 
                K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,surfaceCorrection,T,ndof,par_omega,c,model,damage,history);
            else
                K = analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history);
            end
        end
        disp('Stiffness matrix done.')
        %bn = bn*V;
        du = -K\(bn*V);
        disp("Solution found for the step " + int2str(n) + " out of " + int2str(n_tot))
        un(:,n) = u_trial + du;
        [r_vec,history,phi(:,n)] = getForce(x,un(:,n),T,bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V,par_omega,c,model,damage,history); % Update to include arbitrary displacement kinematic conditions
        r = norm(r_vec,Inf);
    end
    % Energy
    if model.dilatation
        [theta,history.theta] = dilatation(x,un(:,n),family,partialAreas,surfaceCorrection,[],idb,par_omega,c,model,damage,history,0); 
    end
    for ii=1:length(x)
       dofi = dof_vec(ii,:);
       if model.dilatation
        energy.W(ii,n) = strainEnergyDensity(x,un(:,n),theta,family(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,c,model,damage,history.S(ii,:),history.theta)*V;
       else 
        energy.W(ii,n) = strainEnergyDensity(x,un(:,n),[],family(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,c,model,damage,history.S(ii,:))*V;
       end
       if n>1
       bn_1 = b*((n-1)/n_tot); % b_(n-1)
       energy.EW(ii,n) = dot(bn(dofi)+bn_1(dofi),un(dofi,n)-un(dofi,n-1))/2*V + energy.EW(ii,n-1);
       else
       energy.EW(ii,n) = dot(bn(dofi)-0,un(dofi,n)-0)/2*V;
       end
    end
end

end

function [f,history,phi] = getForce(x,u,T,b,familyMat,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V,par_omega,c,model,damage,history)
N = length(x);
f = zeros(size(u));
% Evaluate dilatation
if model.dilatation
    [theta,history.theta] = dilatation(x,u,familyMat,partialAreas,surfaceCorrection,[],idb,par_omega,c,model,damage,history,0);
end
% {Evaluate the force, history variables and damage}
phi = zeros(length(x),1);
for ii = 1:N
    family = familyMat(ii,familyMat(ii,:)~=0); % Neighbours node list
    dofi = dof_vec(ii,:);
    neig_index = 1:length(family);
    jj = family;
    noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
        if model.dilatation
            [fij,history.S(ii,neig_index),mu_j] = T(x,u,theta,ii,jj,dof_vec,par_omega,c,model,[],damage,0,history.S(ii,neig_index),history.theta,noFail);
        else
            [fij,history.S(ii,neig_index),mu_j] = T(x,u,ii,jj,dof_vec,par_omega,c,model,[],damage,[],history.S(ii,neig_index),noFail);
        end
        Vj = partialAreas(ii,neig_index)';
        lambda = surfaceCorrection(ii,neig_index)';
        f_i = sum(fij.*Vj.*lambda);
        f(dofi) = f_i';
        % Damage index
        areaTot = sum(Vj);
        partialDamage = sum(mu_j.*Vj);
        phi(ii) = 1 - partialDamage/areaTot;
end
f = (f + b)*V; % N
% Change it to add boundary conditions
penalty = 1e10;
if ~isempty(bc_set)
    f(ndof+1:end) = penalty*bc_set(:,2); % The second collumn of bc_set represents the value of the constrain
end
end