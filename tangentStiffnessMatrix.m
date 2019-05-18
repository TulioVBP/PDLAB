%% IMPLEMENTATION OF QUASI-STATICS SOLVER - Algorithm

function K = tangentStiffnessMatrix(x,u,idb,family,partialAreas,T,ndof,par_omega,c)
global model
h = norm(x(1,:) - x(2,:)); % Nodal spacing
epsilon = h*1e-6; % According to the roadmap
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
    transvList = [ii family(ii,family(ii,:)~=0)]; % Node i and all neighbours of node i
    for jj= transvList
        %% Evaluate the force state at x1 under perturbations of displacement
        dofj = dof_vec(jj,:);
        for rr = dofj%(2*jj-1):2*jj % For each displacement degree of freedom of node jj: (2*jj-1) = e1 and 2*jj = e2
            epsilon_vector(rr) = epsilon;
            u_plus = u + epsilon_vector;
            u_minus = u - epsilon_vector;
            for kk = family(ii,family(ii,:)~=0)
                dofi = dof_vec(ii,:);
                %[T_plus,~,~] = T(x(ii,:),x(kk,:),u_plus(dofi)',u_plus(dofk)');
                if model.dilatation
                    neigh_Index = find(family(ii,:)==kk);
                    [T_plus,~,~] = T(x,u_plus,ii,dof_vec,family,partialAreas,neigh_Index,par_omega,c);
                    [T_minus,~,~] = T(x,u_minus,ii,dof_vec,family,partialAreas,neigh_Index,par_omega,c);
                else
                    [T_plus,~,~] = T(x,u_plus,ii,kk,dof_vec,par_omega,c);
                    [T_minus,~,~] = T(x,u_minus,ii,kk,dof_vec,par_omega,c);
                end
                %[T_minus,~,~] = T(x(ii,:),x(kk,:),u_minus(dofi)',u_minus(dofk)');
                f_plus = T_plus*partialAreas(ii,family(ii,:)==kk)*h^2; % S_max set to zero
                f_minus = T_minus*partialAreas(ii,family(ii,:)==kk)*h^2; % S_max set to zero again
                f_diff = f_plus - f_minus;
                %for ss = dofk % For each displacement degree of freedom of node kk: (2*jj-1) = e1 and 2*jj = e2
                    K(dofi,rr) = K(dofi,rr) + f_diff'/2/epsilon; % ss+2*(1-kk) returns 1 or 2, depending if ss is the first or the second
                %end
            end
            epsilon_vector(rr) = 0; % Resetting to zero
        end
    end
    disp("Constructing the stiffness matrix: " + num2str(ii/length(x)*100) + "%")
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
