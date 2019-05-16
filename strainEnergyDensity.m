function W = strainEnergyDensity(x,u,family,partialAreas,ii,idb,history)
global c1 c2 horizon omega model
    W = 0;
    familySet = family(family~=0);
    dofi = [idb(2*ii-1) idb(2*ii)];
    m = weightedVolume(horizon,omega);
    % Evaluate dilatation
    if model.dilatation
        theta_i = dilatation(x,u,family,partialAreas,ii,idb,m);
    end
    neigh_ind = 1;
    for jj = familySet
            xi = x(jj,:) - x(ii,:);
            dofj = [idb(2*jj-1) idb(2*jj)];
            eta = u(dofj)' - u(dofi)';
            norma = norm(xi);
            s = (norm(xi+eta) - norma)/norma;
            switch model.name
                case "PMB"
                    if exist('history','var')
                        mu = damageFactor(history(neigh_ind),x(ii,:),x(jj,:)); % NoFail not required
                    else
                        mu = 1;
                    end
                    p = antiderivative(s);
                    w = 1/2*c1*influenceFunction(norma,horizon,omega)*norma^2*p*mu;
                    W = W + w*partialAreas(neigh_ind);
                case "Linearized LPS bond-based"
                    extension = dot(eta,xi)/norma;
                    w = 1/2*c1*influenceFunction(norma,horizon,omega)*extension^2/2;
                    W = W + w*partialAreas(neigh_ind);
                case "Linearized LPS"
                    kappa = c1*m/2+c2*m/3;
                    alfa = c2;
                    w = alfa/2*influenceFunction(norma,horizon,omega)*norma^2*(dot(eta,xi)/norma^2 - theta_i/3)^2;
                    W = W + w*partialAreas(neigh_ind);
                    if jj == familySet(end)
                        W = W + kappa*theta_i^2/2;
                    end
                case "Lipton Free Damage"
                    if exist('history','var')
                        H = damageFactor(history(1,neigh_ind,:),x(ii,:),x(jj,:));
                    else
                        H = [1 1 1];
                    end
                    Slin = dot(xi,eta)/norma^2;
                    V_delta = pi*horizon^2;
                    W = W+2/V_delta*(influenceFunction(norma,horizon,omega)*norma/horizon * H(1)*f_potential(Slin*sqrt(norma)))*partialAreas(neigh_ind);
                    if jj == familySet(end)
                        W = W + 1/horizon^2*H(2)*g_potential(theta_i);
                    end
                case "LPS 2D"
                    elong = norm(xi+eta) - norm(xi);
                    W =  W + c2/2*influenceFunction(norma,horizon,omega)*(elong-theta_i*norma/3)^2*partialAreas(neigh_ind);
                    if jj == familySet(end)
                        W = W + c1*theta_i^2/2;
                    end
                otherwise
                    break
            end
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

function theta = dilatation(x,u,family,partialAreas,ii,idb,m)
    global nu model horizon omega
    familySet = family(family~=0);
    neigh_ind = 1;
    theta = 0;
    dofi = [idb(2*ii-1) idb(2*ii)];
    for jj = familySet
        dofj = [idb(2*jj-1) idb(2*jj)];
        xi = x(jj,:)-x(ii,:);
        norma = norm(xi);
        eta = u(dofj)'-u(dofi)';
        switch model.name            
            case "Linearized LPS"
                theta= theta + 3/m*influenceFunction(norma,horizon,omega)*dot(eta,xi)*partialAreas(neigh_ind);
            case "Lipton Free Damage"
                V_delta = pi*horizon^2;
                S_linear = dot(xi,eta)/norma^2;
                theta = theta + 1/V_delta*influenceFunction(norma,horizon,omega)*norma^2*S_linear*partialAreas(neigh_ind);                
            case "LPS 2D"
                elong = norm(xi+eta) - norm(xi);
                theta = theta + 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,horizon,omega)*norma*elong*partialAreas(neigh_ind);
            otherwise
                break;
        end
        neigh_ind = neigh_ind+1;
    end
end

function ff = f_potential(x)
global damageOn c1
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x <= r1
        ff = c1*x^2/2;
    elseif x > r2
        ff = x;
    end
else
    ff = c1*x^2/2;
end
end

function gg = g_potential(x)
global damageOn c2
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x <= r1
        gg = c2*x^2/2;
    elseif x > r2
        gg = x;
    end
else
    gg = c2*x^2/2;
end

end
