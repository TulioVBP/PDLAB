function [ndof,idb,bc_set,b] = boundaryCondition(x,stresses,m,h,A)
%% INPUT parameter: 
% - x: the nodes position 
%% OUTPUT parameter
% - ndof: Number of degree of freedoms
% - idb: collumn vector that index each of its row (corresponding to the
%        matrix dof) to the actual dof
% EXAMPLE: idb = [1 3 4 6 7 8] means that the second row of the stiffness
%                matrix is related to the 3rd dof freedom of the system
    %% DEFINE THE BOUNDARY CONDITIONS
    bc_set = []; % Each degree of freedom have a index number related
    %% DEFINE THE IDB VECTOR
    ndof = 2*length(x) - length(bc_set);
    %in_ndof = 1:2*length(x);
    idb = zeros(2*length(x),1);
    %ii = 1;
%     for jj = 1:length(in_ndof)
%            check = sum(bc_set == in_ndof(jj));
%            if check ~= 1
%                idb(ii) = in_ndof(jj);
%                ii = ii +1;
%            end
%     end
    id_dof = 1;
    id_const = ndof+1;
    for ii=1:2*length(x)
        check = sum(bc_set == ii);
         if check ~= 1
             % Free degree of freedom
             idb(ii) = id_dof;
             id_dof = id_dof + 1;
         else
             % Constrained degree of freedom
             idb(ii) = id_const;
             id_const = id_const + 1;
         end
    end
    %% DEFINE THE BODY FORCE
    b_old = bodyForce(x,stresses, m, h, A);
    b = zeros(2*length(x),1);
    for ii = 1:length(b_old)
        dofi = [find(idb == (2*ii-1)) find(idb == (2*ii))];
        b(dofi) = b_old(ii,:)'; % Assigning the body force values to the collumn vector
    end
end