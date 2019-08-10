function [f,history,mu] = interactionForce_PMB(x,u,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
%% INPUT
% x - node position matrix
% u - degree of freedom displacement vector
% ii - index of the i-th node
% jj - index of the j-th node inside i's family
% dof_vec - matrix with the degrees of freedom corresponding for each node
% separatorDamage - doesn't do anything but is useful to separate the
%                   normal variables to the ones needed for damage simulation
% dt - step time
% history - maximum stretch for the given bond
% nofail - boolean variable that takes true if 
%% OUTPUT
% f: vector state force between j and i nodes
% history: maximum stretch for each bond
%% CODE
    x_i = x(ii,:); x_j = x(jj,:);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    u = u'; % Nx1 -> 1xN
    u_i = u(dofi); u_j = u(dofj);
    xi = x_j - x_i; % \xi
    eta = u_j - u_i; % \eta
    norma = vecnorm(xi')'; 
    S = (vecnorm(eta'+xi')' - norma)./norma; % Calculate stretch
    ee = (eta + xi)./vecnorm(eta'+xi')'; % Versor
    % Updating maximum stretch
    %if exist('history','var') ~=0
    if nargin > 10  && damage.damageOn% Damage considered
        S_max = history';
        history(S>S_max) = S(S>S_max);
        S_max = history';
        % Evaluating the damage factor
        mu = damageFactor(S_max,ii,1:length(jj),damage,noFail,model); % If noFail is true then we will always have mu as one
        if mu < 1
            a = 1;
        end
    else % No damage considered
        history = S;
        mu = ones(length(S),1);
        noFail = [];
    end
    % Evaluating the force interaction
    f = c(1)*influenceFunction(norma,par_omega).*norma.*fscalar(S.*mu,damage,noFail).*ee; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi|  
end

function ff = fscalar(x,damage,noFail)
%%
%global S0 S1 damageOn
if damage.damageOn
    % Damage dependent crack
    alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
    if damage.phi > alfa
        Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
    else
        Sc = damage.Sc;
    end
    
    % NEW FORMULATION
    ff = (x<Sc).*x;
    ff(noFail) = x(noFail);
else
    ff = x;
end
end
