function theta = dilatation(x,u,family,partialAreas,surfaceCorrection,transvList,idb,par_omega,c,model)
    horizon = par_omega(1);
    m = weightedVolume(par_omega);
    if isempty(transvList) % Not a specific range of nodes was chosen
       transvList = 1:length(x);
    end
    theta = zeros(length(transvList),1); % Initialize dilatation vector
    transv_ind = 1; % Introduced so that we can pass as argument a smaller matrix
    for ii = transvList
        dofi = [idb(2*ii-1) idb(2*ii)];
        familySet = family(transv_ind,family(transv_ind,:)~=0);
        neigh_ind = 1;
        for jj = familySet
            dofj = [idb(2*jj-1) idb(2*jj)];
            xi = x(jj,:)-x(ii,:);
            norma = norm(xi);
            eta = u(dofj)'-u(dofi)';
            switch model.name            
                case "Linearized LPS"
                    theta(transv_ind) = theta(transv_ind) + 3/m*influenceFunction(norma,par_omega)*dot(eta,xi)*partialAreas(transv_ind,neigh_ind)*surfaceCorrection(transv_ind,neigh_ind);
                case "Lipton Free Damage"
                    V_delta = pi*horizon^2;
                    S_linear = dot(xi,eta)/norma^2;
                    theta(transv_ind) = theta(transv_ind) + 1/V_delta*influenceFunction(norma,par_omega)*norma^2*S_linear*partialAreas(transv_ind,neigh_ind)*surfaceCorrection(transv_ind,neigh_ind);                
                case "LPS 2D"
                    nu = c(3);
                    elong = norm(xi+eta) - norm(xi);
                    theta(transv_ind) = theta(transv_ind) + 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,par_omega)*norma*elong*partialAreas(transv_ind,neigh_ind)*surfaceCorrection(transv_ind,neigh_ind);
                otherwise
                    break;
            end
            neigh_ind = neigh_ind+1;
        end
        transv_ind = transv_ind + 1;
    end
end
