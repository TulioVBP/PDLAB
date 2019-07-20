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
            w = 1/V_delta*(influenceFunction(norma,par_omega).*norma/horizon.*H(:,1).*f_potential(Slin,sqrt(norma),c,damage));
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)') + 1/horizon^2*H(1,2)*g_potential(theta_i,c,damage);
        case 4 %"LPS 2D"
            if nargin > 11 && damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                mu = damageFactor(historyS(neigh_ind)',ii,neigh_ind,damage,noFail,model); % NoFail not required
                p = antiderivativePMB(s,damage,noFail);
            else
                mu = ones(length(jj),1);
                p = antiderivativePMB(s,damage,false);
            end
            %elong = vecnorm(xi'+eta')' - norma;
            nu = c(3);
            %w = c(2)/2*influenceFunction(norma,par_omega).*elong.^2;
            w = c(2)/2*influenceFunction(norma,par_omega).*norma.^2.*(2*p).*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)')+ (c(1)/2 + c(2)*m/3*(1/6 - (nu-1)/(2*(2*nu-1))))*theta_i^2;
        case 5 %"PMB"
            if nargin > 11 && damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                mu = damageFactor(historyS(neigh_ind)',ii,neigh_ind,damage,noFail,model); % NoFail not required
                p = antiderivativePMB(s,damage,noFail);
            else
                mu = ones(length(jj),1);
                p = antiderivativePMB(s,damage,false);
            end
            w = 1/2*c(1)*influenceFunction(norma,par_omega).*norma.^2.*p.*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)');
        otherwise
    end
    %neigh_ind = neigh_ind + 1;
    %end
end

function p = antiderivative(x,damage,noFail)
    % Modified PMB model
    if damage.damageOn
        % Damage dependent crack
        alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
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


function p = antiderivativePMB(x,damage,noFail)
    % Modified PMB model
    if damage.damageOn
        % Damage dependent crack
        alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
        if damage.phi > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
        else
            Sc = damage.Sc;
        end
        p = (x<Sc).*x.^2/2;
        % {Correcting the noFail}
        p(noFail) = x(noFail).^2/2; 
    else
        p = x.^2/2;
    end
end

function ff = f_potential(S,norma,c,damage)
r1 = damage.Sc.*norma;
%r1 = 3; % Uncomment for a better result
r2 = r1;
x = S.*norma;
if damage.damageOn
    ff = (x <= r1).*c(1).*x.^2/2 + (x > r2).*x;
else
    ff = c(1)*x.^2/2;
end
end

function gg = g_potential(x,c,damage)
r1 = damage.thetaC;
r2 = r1;
if damage.damageOn
    gg = (x <= r1).*c(2).*x.^2/2 + (x > r2).*x;
else
    gg = c(2)*x^2/2;
end

end
