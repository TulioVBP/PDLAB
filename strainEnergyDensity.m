function W = strainEnergyDensity(x,u,theta,family,partialAreas,surfaceCorrection,ii,idb,par_omega,c,model,damage,historyS,historyT)
    W = 0;
    familySet = family(family~=0);
    dofi = [idb(2*ii-1) idb(2*ii)];
    m = weightedVolume(par_omega);
    horizon = par_omega(1);
    % Evaluate dilatation
    if model.dilatation
        %theta_i = dilatation(x,u,family,partialAreas,ii,idb,m,par_omega,c,model);
        theta_i = theta(ii);
    end
    neigh_ind = 1:length(familySet);
    jj = familySet;
    xi = x(jj,:) - x(ii,:); % Nx2
    dofj = [idb(2*jj-1) idb(2*jj)]; % NX2
    u = u';
    eta = u(dofj) - u(dofi); 
    norma = vecnorm(xi')';
    s = (vecnorm(xi'+eta')' - norma)./norma;
    switch model.number
        case 1 % PMB
            if nargin > 11 && damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                mu = damageFactor(historyS(neigh_ind)',ii,neigh_ind,damage,noFail,model); % NoFail not required
                p = antiderivative(s,damage,noFail);
            else
                mu = ones(length(jj),1);
                p = antiderivative(s,damage,false);
            end
            w = 1/2*c(1)*influenceFunction(norma,par_omega).*norma.^2.*p.*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)');
        case 2 % "Linearized LPS bond-based"
            extension = dot(eta',xi')'./norma;
            w = 1/2*c(1)*influenceFunction(norma,par_omega).*extension.^2/2;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)');
        case 6 %"Linearized LPS"
            kappa = c(1)*m/2+c(2)*m/3;
            alfa = c(2);
            w = alfa/2*influenceFunction(norma,par_omega)*norma^2*(dot(eta,xi)/norma^2 - theta_i/3)^2;
            W = W + w*partialAreas(neigh_ind)*surfaceCorrection(neigh_ind);
            if b_familyEnd
                W = W + kappa*theta_i^2/2;
            end
        case 3 %"Lipton Free Damage"
            if nargin > 11 && damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                XX = [historyS(neigh_ind)', historyT(ii)*ones(length(jj),1),historyT(jj)];
                H = damageFactor(XX,ii,neigh_ind,damage,noFail,model);
            else
                H = ones(length(jj),3);
            end
            Slin = dot(xi',eta')'./norma.^2;
            V_delta = pi*horizon^2;
            w = 1/V_delta*(influenceFunction(norma,par_omega).*norma/horizon.*H(:,1).*f_potential(Slin.*sqrt(norma),c,damage));
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)') + 1/horizon^2*H(1,2)*g_potential(theta_i,c,damage);
        case 4 %"LPS 2D"
            elong = vecnorm(xi'+eta')' - norma;
            nu = c(3);
            %W =  W + c(2)/2*influenceFunction(norma,par_omega)*(elong-theta_i*norma/3)^2*partialAreas(neigh_ind)*surfaceCorrection(neigh_ind);
            w = c(2)/2*influenceFunction(norma,par_omega).*elong.^2;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)')+ (c(1)/2 + c(2)*m/3*(1/6 - (nu-1)/(2*(2*nu-1))))*theta_i^2;
%             if b_familyEnd
%                 %W = W + c(1)*theta_i^2/2;
%                 W = W + (c(1)/2 + c(2)*m/3*(1/6 - (nu-1)/(2*(2*nu-1))))*theta_i^2;
%             end
        otherwise
    end
    %neigh_ind = neigh_ind + 1;
    %end
end

function p = antiderivative(x,damage,noFail)
    % Modified PMB model
    if damage.damageOn
        % Damage dependent crack
        alfa = 0.2; beta = 0.2; gamma = 1.4;
        if damage.phi > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
        else
            Sc = damage.Sc;
        end
        S0 = [-0.98 0.95*Sc]; % S0- and S0+
        S1 = [-0.99 1.05*Sc]; % S1- and S1+
        % Evaluate integration constants
        A = [1 0 0 0 0; 0 1 0 0 0; 1 0 0 -1 0; 0 0 1 0 0; 0 0 -1 0 1];
        b = [S0(1)^2/2 - S0(1)/(S0(1) - S1(1))*(S0(1)^2/2 - S1(1)*S0(1));
            0;
            S0(1)/(S0(1) - S1(1))*(S1(1)^2/2);
            S0(2)^2/2 - S0(2)/(S1(2) - S0(2))*(-S0(2)^2/2 + S1(2)*S0(2));
            S0(2)/(S1(2) - S0(2))*S1(2)^2/2];
        C = A\b;
        p = (x<=S1(1)).* C(4) + (x<=S0(1)).*(x>S1(1)).*(S0(1)/(S0(1) - S1(1)).*(x.^2/2 - S1(1)*x) + C(1)) ...
            + (x<=S0(2)).*(x>S0(1)).*(x.^2/2 + C(2)) + (x<=S1(2)).*(x>S0(2)).*(S0(2)/(S1(2) - S0(2)).*(S1(2)*x - x.^2/2) + C(3)) ...
            + (x>S1(2)).*C(5);
        % {Correcting the noFail}
        p(noFail) = x(noFail).^2/2; 
    else
        p = x.^2/2;
    end
end

function theta = dilatation(x,u,family,partialAreas,ii,idb,m,par_omega,c,model)
    horizon = par_omega(1);
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
                theta= theta + 3/m*influenceFunction(norma,par_omega)*dot(eta,xi)*partialAreas(neigh_ind);
            case "Lipton Free Damage"
                V_delta = pi*horizon^2;
                S_linear = dot(xi,eta)/norma^2;
                theta = theta + 1/V_delta*influenceFunction(norma,par_omega)*norma^2*S_linear*partialAreas(neigh_ind);                
            case "LPS 2D"
                nu = c(3);
                elong = norm(xi+eta) - norm(xi);
                theta = theta + 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,par_omega)*norma*elong*partialAreas(neigh_ind);
            otherwise
                break;
        end
        neigh_ind = neigh_ind+1;
    end
end

function ff = f_potential(x,c,damage)
r1 = 3.0;
r2 = 3.0;
if damage.damageOn
    ff = (x <= r1).*c(1).*x.^2/2 + (x > r2).*x;
else
    ff = c(1)*x.^2/2;
end
end

function gg = g_potential(x,c,damage)
r1 = 3.0;
r2 = 3.0;
if damage.damageOn
    gg = (x <= r1).*c(2).*x.^2/2 + (x > r2).*x;
else
    gg = c(2)*x^2/2;
end

end
