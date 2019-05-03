%QUASI-STATIC SOLVER

function [un,r] = solver_QuasiStatic(x,n_tot,idb,b,bc_set,family,partialAreas,T,ndof,V)
global model
N = length(idb);
% Step 1 - Initialization
un = zeros(N,n_tot); % N= 2*nn
% Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
for n = 1:n_tot
    % Step 2 - Update the load step n <- n + 1 and pseudo-time t. Update the boundary conditions.
    bn = b*(n/n_tot); % Partial load
    % Step 3 - Evaluate the residual vector, r, and residual r. Determine the convergence criterion
    %          for the load step.
    epsilon = 10^-4;
    r_vec = getForce(x,un(:,n),T,bn,family,partialAreas,dof_vec,ndof,bc_set,V); % Update to include arbitrary displacement kinematic conditions
    r_max = epsilon*max(norm(bn*V,Inf),norm(r_vec+bn*V,Inf)); % Normalizing the maximum residual
    % Step 4 - Assign an initial guess to the trial displacement utrial (for example, utrial = un).
    u_trial = un(:,n);
    % Step 5 - Apply Newton's method to minimize the residual.
    r = norm(r_vec,Inf);
    % -------------------- Newton's method ----------------
    if ~model.linearity
        % Suitable for non-linear models
        while r > r_max
            if ~model.stiffnessAnal 
                K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,T,ndof);
            else
                K = analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,V);
            end
            disp('Stiffness matrix done.')
            du = -K\r_vec;
            disp('Incremental solution found.')
            u_trial = u_trial + du;
            r_vec = getForce(x,u_trial,T,bn,family,partialAreas,dof_vec,ndof,bc_set,V); % Update to include arbitrary displacement kinematic conditions
            r = norm(r_vec,Inf);
            r_max = epsilon*max(norm(bn*V,Inf),norm(r_vec+bn*V,Inf));
            disp("Current residual equal to "+ num2str(r) + ". Maximum residual to be " + num2str(r_max))
        end
        disp("Solution found for the step " + int2str(n) + " out of " + int2str(n_tot))
        un(:,n) = u_trial; 
        un(:,n+1) = u_trial; % For the next iteration
    else
        % Linear model
        if ~model.stiffnessAnal 
            K = tangentStiffnessMatrix(x,u_trial,idb,family,partialAreas,T,ndof);
        else
            K = analyticalStiffnessMatrix(x,u_trial,ndof,idb,family,partialAreas,V);
        end        
        disp('Stiffness matrix done.')
        bn = bn*V;
        du = -K\bn;
        disp('Solution found.')
        un(:,n) = u_trial + du;
    end    
end

end

function f = getForce(x,u,T,b,familyMat,partialAreas,dof_vec,ndof,bc_set,V)
N = length(x);
f = zeros(size(u));
for ii = 1:N
    family = familyMat(ii,familyMat(ii,:)~=0); % Neighbours node list
    dofi = dof_vec(ii,:);
    u_i = u(dofi)';
    neighIndex = 1;
    for jj = family
        dofj = dof_vec(jj,:);
        u_j = u(dofj)';
        %[fij,~] = T(x(ii,:),x(jj,:),u_i,u_j,0,[]); % Density vector force
        [fij,~] = T(x,u,ii,dof_vec,familyMat,neighIndex);
        Vj = partialAreas(ii,familyMat(ii,:) == jj);
        f(dofi) = f(dofi) + (fij')*Vj; % N/m^3
        neighIndex = neighIndex + 1;
    end   
end
f = (f + b)*V; % N
% Change it to add boundary conditions
penalty = 1e10;
f(ndof+1:end) = penalty*bc_set(:,2); % The second collumn of bc_set represents the value of the constrain
end