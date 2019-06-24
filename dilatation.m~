function theta = dilatation(x,u,family,partialAreas,XJ,YJ,idb,par_omega,c,model)
    horizon = par_omega(1);
    m = weightedVolume(par_omega);
    theta = zeros(length(x),1);
    for ii = 1:length(x)
    familySet = family(ii,family(ii,:)~=0);
    neigh_ind = 1;
    theta(ii) = 0;
    dofi = [idb(2*ii-1) idb(2*ii)];
    for jj = familySet
        dofj = [idb(2*jj-1) idb(2*jj)];
        x_j = [XJ(ii,neigh_ind) YJ(ii,neigh_ind)];
        xi = x_j-x(ii,:);
        norma = norm(xi);
        eta = u(dofj)'-u(dofi)';
        switch model.name            
            case "Linearized LPS"
                theta(ii) = theta(ii) + 3/m*influenceFunction(norma,par_omega)*dot(eta,xi)*partialAreas(ii,neigh_ind);
            case "Lipton Free Damage"
                V_delta = pi*horizon^2;
                S_linear = dot(xi,eta)/norma^2;
                theta(ii) = theta(ii) + 1/V_delta*influenceFunction(norma,par_omega)*norma^2*S_linear*partialAreas(ii,neigh_ind);                
            case "LPS 2D"
                nu = c(3);
                elong = norm(xi+eta) - norm(xi);
                theta(ii) = theta(ii) + 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,par_omega)*norma*elong*partialAreas(ii,neigh_ind);
            otherwise
                break;
        end
        neigh_ind = neigh_ind+1;
    end
    end
end