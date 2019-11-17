% Function to generate the stiffness matrix for the quasi-static solver
function A = analyticalStiffnessMatrix(x,u,ndof,idb,familySet,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history)
%% INPUT:
% ------------------------------------------------------------
% - x: position of the nodes
% - u: displacement of the nodes
% - ndof: number of degrees of freedom
% - horizon: peridynamic horizon of the simulation
% - familySet: index of every node j (column) inside i (line) node's family
% - partialArea: partial areas of node in j collumn to the ith node
% - V: scalar volume for each node
% - par_omega: horizon omega alfa
% - c: material's constant
%% OUTPUT:
% ------------------------------------------------------------
% - A: PD static matrix
    penalty = 1e10;
    N = size(x,1);
    A = zeros(2*N,2*N); % 2N GDLs
    m = weightedVolume(par_omega)*ones(length(x),1);
    switch model.name
        case "Linearized LPS bond-based"
            c = c(1)/2*weightedVolume(par_omega);
            for ii = 1:N
                %node_i = ceil(ii/2); % Finding the node related to the
                dofi = [idb(2*ii-1) idb(2*ii)];
                family = familySet(ii,familySet(ii,:)~=0);
                iII = 1;
                for jj = family % j sum
                    dofj = [idb(2*jj-1) idb(2*jj)];
                    xi = x(jj,:) - x(ii,:);
                    normaj = norm(xi);
                    omegaj = influenceFunction(normaj,par_omega);
                    % U
                    if dofi(1) <= ndof
                        % First dof of node ii is free
                        ti1u = c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(1))*xi(1)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii); % Aii
                        ti2u = c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(2))*xi(1)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii); % Aip
                        tj1u = -c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(1))*xi(1)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);% Aij
                        tj2u = -c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(2))*xi(1)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);% Aijp
                        A(dofi(1),dofi(1)) = A(dofi(1),dofi(1)) + ti1u;
                        A(dofi(1),dofj(1)) = tj1u;
                        A(dofi(1),dofi(2)) = A(dofi(1),dofi(2)) + ti2u;
                        A(dofi(1),dofj(2)) = tj2u;
                    else
                        % Constraint nodes
                        A(dofi(1),dofi(1)) = penalty;
                    end
                    if dofi(2) <= ndof
                        % V
                        ti1v = c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(1))*xi(2)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);
                        ti2v = c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(2))*xi(2)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);
                        tj1v = -c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(1))*xi(2)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);
                        tj2v = -c*(1/m(jj) + 1/m(ii))*(omegaj/(normaj^2)*xi(2))*xi(2)*partialAreas(ii,iII)*surfaceCorrection(ii,iII)*V(ii);
                        A(dofi(2),dofi(1)) = A(dofi(2),dofi(1)) + ti1v;
                        A(dofi(2),dofj(1)) = tj1v;
                        A(dofi(2),dofi(2)) = A(dofi(2),dofi(2)) + ti2v;
                        A(dofi(2),dofj(2)) = tj2v;
                    else
                        % Constraint nodes
                        A(dofi(2),dofi(2)) = penalty;
                    end

                    % Upload the neigh index
                    iII = iII + 1;
                end
            end
            A = -A;
        otherwise
            disp("Model not implemented.")
    end
end