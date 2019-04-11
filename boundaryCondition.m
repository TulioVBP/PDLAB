function [ndof,idb,bc_set,b,noFail] = boundaryCondition(x,stresses,m,h,A)
%% INPUT parameter: 
% - x: the nodes position 
% - stresses: [sigma_x, sigma_y, tau_xy]
% - m: number of nodes inside delta
% - h: grid spacing
% - A: element area
%% OUTPUT parameter
% - ndof: Number of degree of freedoms
% - idb: collumn vector that index each of its row (corresponding to the
%        actual dof) to the matrix dof
%        EXAMPLE: idb = [1 3 4 6 7 8] means that the thirf row of the stiffness
%                matrix is related to the 2nd dof freedom of the system
% - bc_set: contains the set of constrained degrees of freedom on the first
%           collumn and its corresponding value on the second collumn
% - b: the actual body force vector
% - noFail: set of nodes for which we have no fail condition (mu = 1 always)
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
    [b_old,noFail] = bodyForce(x,stresses, m, h, A);
    b = zeros(2*length(x),1);
    for ii = 1:length(b_old)
        dofi = [idb(2*ii-1) idb == (2*ii)];
        b(dofi) = b_old(ii,:)'; % Assigning the body force values to the collumn vector
    end
end