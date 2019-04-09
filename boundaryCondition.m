function [ndof,idb,bc_set] = boundaryCondition(x)
%% INPUT parameter: 
% - x: the nodes position 
%% OUTPUT parameter
% - ndof: Number of degree of freedoms
% - idb: collumn vector that index each of its row (corresponding to the
%        matrix dof) to the actual dof
% EXAMPLE: idb = [1 3 4 6 7 8] means that the second row of the stiffness
%                matrix is related to the 3rd dof freedom of the system
    %% DEFINE THE BOUNDARY CONDITIONS
    bc_set = [1 2 5 7]; % Each degree of freedom have a index number related
    %% DEFINE THE IDB VECTOR
    ndof = 2*length(x) - length(bc_set);
    in_ndof = 1:2*length(x);
    idb = zeros(ndof,1);
    ii = 1;
    for jj = 1:length(in_ndof)
           check = sum(bc_set == in_ndof(jj));
           if check ~= 1
               idb(ii) = in_ndof(jj);
               ii = ii +1;
           end
    end
end