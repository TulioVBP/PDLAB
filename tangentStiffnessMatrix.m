%% IMPLEMENTATION OF QUASI-STATICS SOLVER - Algorithm

function K = tangentStiffnessMatrix(x,u,b,idb,family,partialAreas,T,ngdl)
h = norm(x(1,:) - x(2,:)); % Nodal spacing
epsilon = h*1e-6; % According to the roadmap
epsilon_vector = zeros(ngdl);
% Initialize the tangent stiffness matrix to 0
N = length(x);
%% Initialize the tangent stiffness matrix to zero
K = zeros(2*N,2*N);
%% Transverse each node in the discretization
for ii=1:length(x)
    transvList = [ii family(ii,family(ii,:)~=0)]; % Node i and all neighbours of node i
    for jj= transvList
        %% Evaluate the force state at x1 under perturbations of displacement
        for rr = (2*jj-1):2*jj % For each displacement degree of freedom of node jj: (2*jj-1) = e1 and 2*jj = e2
            epsilon_vector(rr) = epsilon;
            u_plus = u + epsilon_vector;
            u_minus = u - epsilon_vector;
            for kk = family(ii,family(ii,:)~=0)
                f_plus = T(x(ii,:),x(kk,:),u_plus((2*ii-1):2*ii),u_plus((2*kk-1):2*kk),0)*partialArea(ii,family(ii,:)==kk)*h^2; % S_max set to zero
                f_minus = T(x(ii,:),x(kk,:),u_minus((2*ii-1):2*ii),u_minus((2*kk-1):2*kk),0)*partialArea(ii,family(ii,:)==kk)*h^2; % S_max set to zero again
                f_diff = f_plus - f_minus;
                for ss = (2*kk-1):2*kk % For each displacement degree of freedom of node kk: (2*jj-1) = e1 and 2*jj = e2
                    K(ss,rr) = K(ss,rr) + f_diff(ss+2*(1-kk))/2/epsilon; % ss+2*(1-kk) returns 1 or 2, depending if ss is the first or the second
                end
            end
            epsilon_vector(rr) = 0; % Resetting to zero
        end
    end
end
