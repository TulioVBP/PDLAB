%% IMPLEMENTATION OF QUASI-STATICS SOLVER - Algorithm

function K = tangentStiffnessMatrix(x,u,idb,family,partialAreas,surfaceCorrection,T,ndof,par_omega,c,model,damage)
h = norm(x(1,:) - x(2,:)); % Nodal spacing
epsilon = h*1e-7; % According to the roadmap
N = length(u);
epsilon_vector = zeros(N,1);
 % Defining the node's degree of freedom index
    dof_vec = zeros(size(x));
    for kk = 1:length(x)
        dof_vec(kk,:) = [idb(2*kk-1) idb(2*kk)];
    end
%% Initialize the tangent stiffness matrix to zero
K = zeros(N,N);
%% Transverse each node in the discretization
for ii=1:length(x)
    if model.dilatation % The system has dilatation
        transvList = [ii family(ii,family(ii,:)~=0)];
        li = length(transvList);
        neighIndex = 1;
        for jj = family(ii,family(ii,:)~=0)
            for ll = family(jj,family(jj,:)~=0)
                if ~sum(transvList == ll)
                    transvList(li+neighIndex) = ll;
                    neighIndex = neighIndex + 1;
                end
            end
        end
    else % No dilatation
        transvList = [ii family(ii,family(ii,:)~=0)]; % Node i and all neighbours of node i
    end
    for jj= transvList
        %% Evaluate the force state at x1 under perturbations of displacement
        dofj = dof_vec(jj,:);
        for rr = dofj%(2*jj-1):2*jj % For each displacement degree of freedom of node jj: (2*jj-1) = e1 and 2*jj = e2
            epsilon_vector(rr) = epsilon;
            u_plus = u + epsilon_vector;
            u_minus = u - epsilon_vector;
            
            % {Evaluate related dilatation}
            theta_plus = zeros(length(x),1); theta_minus = zeros(length(x),1);
            if model.dilatation
                transvListII = [ii family(ii,family(ii,:)~=0)]; % Transversal list of affected dilatations
                theta_plus(transvListII) = dilatation(x,u_plus,family(transvListII,:),partialAreas(transvListII,:),surfaceCorrection(transvListII,:),transvListII,idb,par_omega,c,model);
                theta_minus(transvListII) = dilatation(x,u_minus,family(transvListII,:),partialAreas(transvListII,:),surfaceCorrection(transvListII,:),transvListII,idb,par_omega,c,model);
            end
            % -- 
            
            for kk = family(ii,family(ii,:)~=0)
                dofi = dof_vec(ii,:);
                %[T_plus,~,~] = T(x(ii,:),x(kk,:),u_plus(dofi)',u_plus(dofk)');
                if model.dilatation
                    %neigh_Index = find(family(ii,:)==kk);
                    [T_plus,~,~] = T(x,u_plus,theta_plus,ii,kk,dof_vec,par_omega,c,model,[],damage);
                    [T_minus,~,~] = T(x,u_minus,theta_minus,ii,kk,dof_vec,par_omega,c,model,[],damage);
%                     [T_plus,~,~] = T(x,u_plus,ii,dof_vec,family,partialAreas,neigh_Index,par_omega,c,model,[],damage);
%                     [T_minus,~,~] = T(x,u_minus,ii,dof_vec,family,partialAreas,neigh_Index,par_omega,c,model,[],damage);
                else
                    [T_plus,~,~] = T(x,u_plus,ii,kk,dof_vec,par_omega,c,model,[],damage);
                    [T_minus,~,~] = T(x,u_minus,ii,kk,dof_vec,par_omega,c,model,[],damage);
                end
                %[T_minus,~,~] = T(x(ii,:),x(kk,:),u_minus(dofi)',u_minus(dofk)');
                f_plus = T_plus*partialAreas(ii,family(ii,:)==kk)*surfaceCorrection(ii,family(ii,:)==kk)*h^2; % S_max set to zero
                f_minus = T_minus*partialAreas(ii,family(ii,:)==kk)*surfaceCorrection(ii,family(ii,:)==kk)*h^2; % S_max set to zero again
                f_diff = f_plus - f_minus;
                %for ss = dofk % For each displacement degree of freedom of node kk: (2*jj-1) = e1 and 2*jj = e2
                    K(dofi,rr) = K(dofi,rr) + f_diff'/2/epsilon; % ss+2*(1-kk) returns 1 or 2, depending if ss is the first or the second
                %end
            end
            epsilon_vector(rr) = 0; % Resetting to zero
        end
    end
    if ii == floor(length(x)/2) || ii == floor(length(x)/4) || ii == 3*floor(length(x)/4)
        disp("Constructing the stiffness matrix: " + num2str(ii/length(x)*100) + "%")
    end
end
%% Adapting for the constrained dofs
for ii = ndof+1:length(u)
    K(ii,:) = zeros(size(K(ii,:)));
    K(ii,ii) = 1e10; % Penalty
end
end
% if sum(dofj > ngdl) == 2
%             % Both the dof are constrained
%             dofj = [];
%         elseif sum(dofj>ngdl) == 1
%             % One of the dof is constrained
%             dofj = dofj(dofj<= ngdl);
%         end
