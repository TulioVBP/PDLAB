function [f,historyS,mu] = interactionForce_StateBased(x,u,theta,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,historyS,historyTheta,noFail)
    nu = c(3);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x(jj,:) - x(ii,:); % \xi
    u = u';
    u_i = u(dofi); u_j = u(dofj);
    eta = u_j - u_i; % \eta
    norma = vecnorm(xi')'; 
    elong = vecnorm(xi'+eta')' - norma; % Calculate elongation - linear
    S = elong./norma;
    ee = (xi+eta)./vecnorm(xi'+eta')'; % Versor
    % Defining non-existing parameters
    if nargin > 10  && damage.damageOn% Damage considered
        S_max = historyS';
        historyS(S>S_max) = S(S>S_max);
        S_max = historyS';
        % Evaluating the damage factor
        mu = damageFactor(S_max,ii,1:length(jj),damage,noFail,model); % If noFail is true then we will always have mu as one
    else % No damage considered
        historyS = S;
        mu = ones(length(S),1);
        noFail = [];
    end
    if ~exist('dt','var')
       dt = 0; 
    end
    % ---- Evaluatin the force ----
    % Dilatation term
    %theta = [theta(ii) theta(jj)];
    m = weightedVolume(par_omega);
    ff = fscalar(S,damage,noFail,ii);
    T_ij = 2*(2*nu-1)/(nu-1)*((c(1) + c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega).*norma/m)*theta(ii) ...
        + c(2)*influenceFunction(norma,par_omega).*norma.*mu.*ff;%.*(elong);
    T_ji = 2*(2*nu-1)/(nu-1)*((c(1) + c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega).*norma/m).*theta(jj) ...
        + c(2)*influenceFunction(norma,par_omega).*norma.*mu.*ff;
    f = (T_ij + T_ji).*ee;
    %mu = 1; % No damage in this model
    %historyS = 0;
end

function ff = fscalar(x,damage,noFail,ii)
%%
%global S0 S1 damageOn
if damage.damageOn
    % Damage dependent crack
    alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
    if damage.phi(ii) > alfa
        Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
    else
        Sc = damage.Sc;
    end
    S0 = [-0.98 0.95*Sc]; % S0- and S0+
    S1 = [-0.99 1.05*Sc]; % S1- and S1+

    % NEW FORMULATION
    ff = (x>S1(1)).*(x<S0(1)).*(S0(1)*(x-S1(1))./(S0(1) - S1(1))) + (x>=S0(1)).*(x<=S0(2)).*x ...
        + (x>S0(2)).*(x<S1(2)).*(S0(2)*(S1(2)-x)./(S1(2)-S0(2)));
    % Adding noFail condition
    ff(noFail) = x(noFail);
else
    ff = x;
end
end