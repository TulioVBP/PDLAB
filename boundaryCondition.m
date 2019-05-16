function [ndof,idb,bc_set,bb,noFail] = boundaryCondition(x,stresses,m,h,A)
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
%           collumn and its corresponding value on the second collumn; the
%           third collumn contains the corresponding dof velocity
% - b: the actual body force vector
% - noFail: set of nodes for which we have no fail condition (mu = 1 always)

    %% DEFINE THE BOUNDARY CONDITIONS
    if true
        % Example of BCs
        b = max(x(:,1));
        rightLay = find(x(:,1) == b);
        bottomLay = find(x(:,2) == 0);
        dof_Right = rightLay*2 - 1; % Constrain the x dof
        dof_Bot = bottomLay*2; % Constrain the y dof
        bc_set(:,1) = [dof_Right; dof_Bot]; % Each degree of freedom have a index number related
        bc_set(:,1) = sort(bc_set(:,1));
        bc_set(:,2:3) = zeros(length(bc_set(:,1)),2);
    else
        a = max(x(:,2));
        bottomLay = find(x(:,2) == 0);
        topLay = find(x(:,2) == a);
        bc_set = [];
    end
    %% DEFINE THE IDB VECTOR
    ndof = 2*length(x) - length(bc_set);
    %in_ndof = 1:2*length(x);
    idb = zeros(2*length(x),1);
    id_dof = 1;
    id_const = ndof+1;
    for ii=1:2*length(x)
        if isempty(bc_set)
            check = 0;
        else
            check = sum(bc_set(:,1) == ii);
        end
         if check == 0
             % Free degree of freedom
             idb(ii) = id_dof;
             id_dof = id_dof + 1;
         else
             % Constrained degree of freedom
             idb(ii) = id_const;
             id_const = id_const + 1;
         end
    end
    % Ensure compatibility between idb for constrained nodes and bc_set -
    % TO BE DONE
    %% DEFINE THE BODY FORCE
    [b_old,noFail] = bodyForce(x,stresses, m, h, A);
    bb = zeros(2*length(x),1);
    for ii = 1:length(b_old)
        dofi = [idb(2*ii-1) idb(2*ii)];
        bb(dofi) = b_old(ii,:)'; % Assigning the body force values to the collumn vector
    end
end