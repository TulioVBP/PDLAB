function W = strainEnergyDensity(x,u,family,partialAreas,ii,S_max,idb)
global c1 horizon omega crackIn model
    W = 0;
    familySet = family(family~=0);
    neigh_ind = 1;
    dofi = [idb(2*ii-1) idb(2*ii)];
    for jj = familySet
            xi = x(jj,:) - x(ii,:);
            dofj = [idb(2*jj-1) idb(2*jj)];
            eta = u(dofj)' - u(dofi)';
            norma = norm(xi);
            s = (norm(xi+eta) - norma)/norma;
            switch model
                case "PMB"
                    mu = damageFactor(S_max(neigh_ind),crackIn,x(ii,:),x(jj,:)); % NoFail not required
                    w = c1*influenceFuntion(norma,horizon,omega)*norma^2*s^2*mu/2;
                otherwise
                    break
            end
            W = W+ w*partialAreas(neigh_ind);
            neigh_ind = neigh_ind + 1;
    end
end