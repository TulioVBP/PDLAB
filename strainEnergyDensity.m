function W = strainEnergyDensity(x,u,family,partialAreas,ii,S_max,idb)
global c1 horizon omega model
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
            switch model.name
                case "PMB"
                    mu = damageFactor(S_max(neigh_ind),x(ii,:),x(jj,:)); % NoFail not required
                    p = antiderivative(s);
                    w = c1*influenceFunction(norma,horizon,omega)*norma^2*p*mu;
                case "Linearized LPS bond-based"
                    extension = dot(eta,xi)/norma;
                    w = c1*influenceFunction(norma,horizon,omega)*extension^2/2;
                otherwise
                    break
            end
            W = W+ 1/2*w*partialAreas(neigh_ind);
            neigh_ind = neigh_ind + 1;
    end
end

function p = antiderivative(x)
    % Modified PMB model
    global S0 S1
    % Evaluate integration constants
    A = [1 0 0 0 0; 0 1 0 0 0; 1 0 0 -1 0; 0 0 1 0 0; 0 0 -1 0 1];
    b = [S0(1)^2/2 - S0(1)/(S0(1) - S1(1))*(S0(1)^2/2 - S1(1)*S0(1));
        0;
        S0(1)/(S0(1) - S1(1))*(S1(1)^2/2);
        S0(2)^2/2 - S0(2)/(S1(2) - S0(2))*(-S0(2)^2/2 + S1(2)*S0(2));
        S0(2)/(S1(2) - S0(2))*S1(2)^2/2];
    C = A\b;
    if x <  S1(1)
        p = C(4); % C4
    elseif x < S0(1)
        p = S0(1)/(S0(1) - S1(1))*(x^2/2 - S1(1)*x) + C(1);
    elseif x < S0(2)
        p = x^2/2 + C(2);
    elseif x < S1(2)
        p = S0(2)/(S1(2) - S0(2))*(S1(2)*x - x^2/2) + C(3);
    else
        p = C(5);
    end
end