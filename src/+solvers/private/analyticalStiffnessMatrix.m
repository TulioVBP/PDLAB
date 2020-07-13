% Function to generate the stiffness matrix for the quasi-static solver
function A = analyticalStiffnessMatrix(x,u,ndof,idb,familySet,partialAreas,surfaceCorrection,V,par_omega,c,model,damage,history,mu)
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

if nargin<14
    mu = cell(length(x),1);
    mu(:) = {1};
elseif isempty(mu)
    mu = cell(length(x),1);
    mu(:) = {1};
end
    penalty = 1e10;
    N = size(x,1);
    A = zeros(2*N,2*N); % 2N GDLs
    m = weightedVolume(par_omega)*ones(length(x),1);
    switch model.name
        case "LBB"
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
        case "PMB" % Linearized
            c = c(1)/2*weightedVolume(par_omega);
            u = u';
            for ii = 1:N
                %node_i = ceil(ii/2); % Finding the node related to the
                dofi = [idb(2*ii-1) idb(2*ii)];
                family = familySet(ii,familySet(ii,:)~=0);
                iII = 1:length(family);
                jj = family'; % j sum
                dofj = [idb(2*jj-1) idb(2*jj)];
                eta = u(dofj) - u(dofi);
                xi = x(jj,:) - x(ii,:);
                M = eta+xi;
                normaj = vecnorm(xi')'; 
                norma_eta = vecnorm(M')';
                omegaj = influenceFunction(normaj,par_omega);
                muj = mu{ii};
                % U
                if dofi(1) <= ndof
                    % First dof of node ii is free
                    ti1u = c*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aii
                    ti2u = c*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aip
                    tj1u = -c*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aij
                    tj2u = -c*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aijp
                    A(dofi(1),dofi(1)) = sum(ti1u);
                    A(dofi(1),dofj(:,1)) = tj1u';
                    A(dofi(1),dofi(2)) = sum(ti2u);
                    A(dofi(1),dofj(:,2)) = tj2u';
               else
                    % Constraint nodes
                    A(dofi(1),dofi(1)) = penalty;
               end
               if dofi(2) <= ndof
                  % V
                   ti1v = c*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   ti2v = c*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   tj1v = -c*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   tj2v = -c*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   A(dofi(2),dofi(1)) = sum(ti1v);
                   A(dofi(2),dofj(:,1)) = tj1v';
                   A(dofi(2),dofi(2)) =sum(ti2v);
                   A(dofi(2),dofj(:,2)) = tj2v';
               else
                    % Constraint nodes
                    A(dofi(2),dofi(2)) = penalty;
               end
            end
            A = -A;
        case "LPS-T"
            u = u';
            nu = c(3);
            AA = 2*(2*nu-1)/(nu-1);
            for ii = 1:N
                c1_hat = (c(1) + c(2)/9*m(ii)*(-nu+2)/(2*nu-1));
                dofi = [idb(2*ii-1) idb(2*ii)];
                family = familySet(ii,familySet(ii,:)~=0);
                Nj = length(family);
                iII = 1:Nj;
                jj = family'; % j sum
                dofj = [idb(2*jj-1) idb(2*jj)];
                eta = u(dofj) - u(dofi);
                xi = x(jj,:) - x(ii,:);
                M = eta+xi;
                normaj = vecnorm(xi,2,2); 
                norma_eta = vecnorm(M,2,2);
                omegaj = influenceFunction(normaj,par_omega);
                muj = mu{ii};
                PSI_ij = M./norma_eta;
                % Parameters
                Vij = partialAreas(ii,iII)';
                c1 = c(2)*omegaj;
                c2 = c1_hat*AA/m(ii).*omegaj.*normaj;
                g = omegaj.*normaj;
                % U and V
                if dofi(1) <= ndof || dofi(2) <= ndof
                    % First dof of node ii is free
  
                    ti1 = sum(-2*c1.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij) + sum(-c2.*AA/m(ii).*g.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj).*sum(Vij.*PSI_ij); % Aii
                    ti2 = sum(-2*c1.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij) + sum(-c2.*AA/m(ii).*g.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj).*sum(Vij.*PSI_ij); % Aip
                    tj1 = 2*c1.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij + c2.*AA/m(ii).*(g.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj)*sum(Vij.*PSI_ij);% Aij
                    tj2 = 2*c1.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij + c2.*AA/m(ii).*(g.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj)*sum(Vij.*PSI_ij);% Aijp
                    
                    for Ij = 1:length(jj)
                        j = jj(Ij);
                        kk  = familySet(j,familySet(j,:)~=0)';
                        iIII = 1:length(kk);
                        dofk = [idb(2*kk-1) idb(2*kk)];
                        eta_k = u(dofk) - u(dofj(Ij,:));
                        xi_k = x(kk,:) - x(j,:);
                        M_k = eta_k+xi_k;
                        normak = vecnorm(xi_k,2,2); 
                        
                        norma_etak = vecnorm(M_k,2,2);
                        omegak = influenceFunction(normak,par_omega);
                        muk = mu{j};
                        PSI_jk = M_k./norma_etak;
                        Vjk = partialAreas(j,iIII)';
                        % Parameters
                        c2k = c1_hat*AA/m(ii).*omegak.*normak;
                        gk = omegak.*normak;
                        
                        tj1(Ij,:) = tj1(Ij,:) - c2(Ij).*sum(AA/m(j).*gk.*PSI_jk(:,1).*Vjk.*surfaceCorrection(j,iIII)'.*muk).*Vij(Ij).*PSI_ij(Ij,:); 
                        tj2(Ij,:) = tj2(Ij,:) - c2(Ij).*sum(AA/m(j).*gk.*PSI_jk(:,2).*Vjk.*surfaceCorrection(j,iIII)'.*muk).*Vij(Ij).*PSI_ij(Ij,:);
                        tk1 = c2(Ij)*AA/m(j).*(gk.*PSI_jk(:,1).*Vjk.*surfaceCorrection(j,iIII)'.*muk)*Vij(Ij).*PSI_ij(Ij,:);
                        tk2 = c2(Ij)*AA/m(j).*(gk.*PSI_jk(:,2).*Vjk.*surfaceCorrection(j,iIII)'.*muk)*Vij(Ij).*PSI_ij(Ij,:);
                        
                        if dofi(1) <= ndof
                            A(dofi(1),dofk(:,1)) = A(dofi(1),dofk(:,1)) + tk1(:,1)' * V(ii);
                            A(dofi(1),dofk(:,2)) = A(dofi(1),dofk(:,2)) + tk2(:,1)' * V(ii);
                        end
                        if dofi(2) <= ndof
                            A(dofi(2),dofk(:,1)) = A(dofi(2),dofk(:,1)) + tk1(:,2)' * V(ii);
                            A(dofi(2),dofk(:,2)) = A(dofi(2),dofk(:,2)) + tk2(:,2)' * V(ii);
                        end 
                    end
                    % U
                    if dofi(1) <= ndof
                        A(dofi(1),dofi(1)) = A(dofi(1),dofi(1)) + ti1(1) *V(ii) ;
                        A(dofi(1),dofj(:,1)) = A(dofi(1),dofj(:,1)) + tj1(:,1)' * V(ii);
                        A(dofi(1),dofi(2)) = A(dofi(1),dofi(2)) + ti2(1) * V(ii);
                        A(dofi(1),dofj(:,2)) = A(dofi(1),dofj(:,2)) + tj2(:,1)'* V(ii);
                    else
                        % Constraint nodes
                        A(dofi(1),dofi(1)) = -penalty;
                    end
                   % V
                   if dofi(2) <= ndof
                        A(dofi(2),dofi(1)) = A(dofi(2),dofi(1)) + ti1(2) *V(ii) ;
                        A(dofi(2),dofj(:,1)) = A(dofi(2),dofj(:,1)) + tj1(:,2)' * V(ii);
                        A(dofi(2),dofi(2)) = A(dofi(2),dofi(2)) + ti2(2) * V(ii);
                        A(dofi(2),dofj(:,2)) = A(dofi(2),dofj(:,2)) + tj2(:,2)' * V(ii);
                   else
                        % Constraint nodes
                        A(dofi(2),dofi(2)) = -penalty;
                   end
                else
                    A(dofi(1),dofi(1)) = -penalty;
                    A(dofi(2),dofi(2)) = -penalty;
                end
            end
        otherwise
            error("Model not implemented.")
    end
end