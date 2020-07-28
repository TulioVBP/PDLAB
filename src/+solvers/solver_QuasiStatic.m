% QUASI-STATIC SOLVER

function [un,r,energy,phi] = solver_QuasiStatic(x,n_tot,idb,b,bc_set,family,partialAreas,surfaceCorrection,model,par_omega,ndof,V,damage,noFailZone)
N = length(idb);
if length(V) == 1
   V = V*ones(length(x),1); 
end
% Create v_dof
V_DOF = zeros(2*length(x),1);
for ii = 1:length(x)
   V_DOF([2*ii-1 2*ii]) = V(ii);  
end
penalty = 1e10;
%% Damage variables 
%{Defining cracking trespassing matrix}
    %if damage.damageOn
        crackSegments = size(damage.crackIn,1); % At least 2
        damage.checkCrack = zeros(size(model.history.S));
        for ii = 1:size(model.history.S)
            x_j = x(family(ii,family(ii,:)~=0),:);
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
    %end
    % {No fail to damage variable}
    damage.noFail = noFailZone;   
    phi = zeros(length(x),n_tot); % Damage index
    

%% Step 1 - Initialization
un = initialU0(N, n_tot,bc_set); % N= 2*nn
energy.W = zeros(length(x),n_tot);
energy.KE = zeros(length(x),n_tot);
energy.EW = zeros(length(x),n_tot);
f_int = zeros(N,n_tot);
history.S = model.history.S;
if model.b_dilatation
    history.theta = model.history.theta;
end
% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end

for n = 1:n_tot
    %% Step 2 - Update the load step n <- n + 1 and pseudo-time t. Update the boundary conditions.
    bn = b*(n/n_tot); % Partial load
    if ~isempty(bc_set)
        bc_setn = [bc_set(:,1),bc_set(:,2)*n/n_tot];
    end
    %% Step 2.5 - Assign an initial guess to the trial displacement utrial (for example, utrial = un).
    u_trial = un(:,n);
    if ~isempty(bc_set)
        u_trial(ndof+1:end) = bc_setn(:,2); % minus is because it is on the right hand side
    end
    %% Step 3 - Evaluate the residual vector, r, and residual r. Determine the convergence criterion
    %          for the load step.
    epsilon = 10^-4;
    damage.phi = phi(:,n); % Accessing current damage situation
    [r_vec,history,phi(:,n),f_int(:,n)] = getForce(x,u_trial,bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_setn,V_DOF,par_omega,model,damage,history); % Update to include arbitrary displacement kinematic conditions
    r_max = epsilon*max(norm(bn(1:ndof).*V_DOF(1:ndof),Inf),norm(f_int(1:ndof,n).*V_DOF(1:ndof),Inf)); % Normalizing the maximum residual
    if r_max == 0 % No forces on the non-constrain nodes
        r_max = 10^-4;
    end
    %% Step 5 - Apply Newton's method to minimize the residual.
    r = norm(r_vec(1:ndof),Inf);
    alpha = 1;
    % -------------------- Newton's method ----------------
    if ~model.b_linearity
        % Suitable for non-linear models
        iter = 1;
        while r > r_max
            % {Damage}
            damage.phi = phi(:,n); % Accessing current damage situation
            if ~model.b_stiffnessAnal 
                K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,V,surfaceCorrection,ndof,par_omega,model,damage,history);
            else
                K = model.analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,damage,history);
            end
            disp('Stiffness matrix done.')
            du = -K\r_vec;
            disp('Incremental solution found.')
            u_trial = u_trial + alpha*du;
            [r_vec,history,phi(:,n),f_int(:,n)] = getForce(x,u_trial,bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_setn,V_DOF,par_omega,model,damage,history); % Update to include arbitrary displacement kinematic conditions
            r = norm(r_vec(1:ndof),Inf);
            %r_max = epsilon*max(norm(bn*V,Inf),norm(r_vec-bn*V,Inf));
            disp("Iter:" + int2str(iter) + " Residual equal to "+ num2str(r) + ". Maximum residual to be " + num2str(r_max))
            iter = iter + 1;
        end
        disp("Solution found for the step " + int2str(n) + " out of " + int2str(n_tot))
        un(:,n) = u_trial;
        if n < n_tot
            un(:,n+1) = u_trial; % For the next iteration
        end
    else
        u_trial = un(:,n);
        % Linear model
        % {Damage}
        damage.phi = phi(:,n); % Accessing current damage situation
        if n < 2 % If the model is linear, there is no need to find the matrix more than once
            if ~model.b_stiffnessAnal 
                K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,V,surfaceCorrection,ndof,par_omega,model,damage,history);
            else
                K = model.analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,surfaceCorrection,V,par_omega,damage,history);
            end
        end
        disp('Stiffness matrix done.')
        %bn = bn*V;
        ff = bn.*V_DOF;
        if ~isempty(bc_set)
            ff(ndof+1:end) = -penalty*bc_set(:,2)*(n/n_tot);
        end
        du = -K\(ff);
        disp("Solution found for the step " + int2str(n) + " out of " + int2str(n_tot))
        un(:,n) = u_trial + du;
        [r_vec,history,phi(:,n),f_int(:,n)] = getForce(x,un(:,n),bn,family,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V_DOF,par_omega,model,damage,history); % Update to include arbitrary displacement kinematic conditions
        r = norm(r_vec(1:ndof),Inf);
    end
    % Energy
    if model.b_dilatation
        [theta,history.theta] = model.dilatation(x,un(:,n),family,partialAreas,surfaceCorrection,[],idb,par_omega,damage,history.S,history.theta); 
    end
    if ~isempty(bc_set)
       con_dof = idb(bc_set(:,1));
    else
       con_dof = [];
    end
    for ii=1:length(x)
       dofi = dof_vec(ii,:);
       if model.b_dilatation
        energy.W(ii,n) = model.strainEnergyDensity(x,un(:,n),theta,family(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,damage,history.S(ii,:),history.theta)*V(ii);
       else 
        energy.W(ii,n) = model.strainEnergyDensity(x,un(:,n),family(ii,:),partialAreas(ii,:),surfaceCorrection(ii,:),ii,idb,par_omega,damage,history.S(ii,:))*V(ii);
       end
       if n>1
           bn_1 = b*((n-1)/n_tot); % b_(n-1)
           % {External work realized by the displacement constraint}
            if sum(dofi(1)==con_dof)
                bn(dofi(1)) = -f_int(dofi(1),n);
                bn_1(dofi(1)) = -f_int(dofi(1),n-1);
                %du = u(con_dof,n) - u(con_dof,n-1);
                %add_ext = bf.*du;
            end
            if sum(dofi(2)==con_dof)
                bn(dofi(2)) = -f_int(dofi(2),n);
                bn_1(dofi(2)) = -f_int(dofi(2),n-1);
            end
           energy.EW(ii,n) = dot(bn(dofi)+bn_1(dofi),un(dofi,n)-un(dofi,n-1))/2*V(ii) + energy.EW(ii,n-1);
       else
           if sum(dofi(1)==con_dof)
                bn(dofi(1)) = -f_int(dofi(1),n);
            end
            if sum(dofi(2)==con_dof)
                bn(dofi(2)) = -f_int(dofi(2),n);
            end
           energy.EW(ii,n) = dot(bn(dofi)-0,un(dofi,n)-0)/2*V(ii);
       end
    end
end

end

function [f,history,phi,f_model] = getForce(x,u,b,familyMat,partialAreas,surfaceCorrection,dof_vec,idb,ndof,bc_set,V,par_omega,model,damage,history)
N = length(x);
f = zeros(size(u));
% Evaluate dilatation
if model.b_dilatation
    [theta,history.theta] = model.dilatation(x,u,familyMat,partialAreas,surfaceCorrection,[],idb,par_omega,damage,history.S,history.theta);
end
% {Evaluate the force, history variables and damage}
phi = zeros(length(x),1);
for ii = 1:N
    family = familyMat(ii,familyMat(ii,:)~=0); % Neighbours node list
    dofi = dof_vec(ii,:);
    neig_index = 1:length(family);
    jj = family;
    noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
        if model.b_dilatation
            [fij,history.S(ii,neig_index),mu_j] = model.T(x,u,theta,ii,jj,dof_vec,par_omega,[],damage,history.S(ii,neig_index),history.theta,noFail);
        else
            [fij,history.S(ii,neig_index),mu_j] = model.T(x,u,ii,jj,dof_vec,par_omega,[],damage,history.S(ii,neig_index),noFail);
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
f_model = f; % Internal force only for all points (including b)
f = (f + b).*V; % 2N
% Change it to add boundary conditions
penalty = 1e10;
if ~isempty(bc_set)
    f(ndof+1:end) = - penalty*zeros(size(bc_set(:,2))); % The second collumn of bc_set represents the value of the constrain; minus for the relation u = -K\f;
end
end

function un = initialU0(N, n_tot,bc_set)
n_con = size(bc_set,1);
un = zeros(N,n_tot);
%un(end-n_con+1:end,1) = bc_set(:,2)/n_tot;
end