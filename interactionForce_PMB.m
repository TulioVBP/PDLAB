function [f,history,mu] = interactionForce_PMB(x,u,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
%% INPUT
% x - node position matrix (changed with quadrature points)
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
    u_i = u(dofi)'; u_j = u(dofj)';
    xi = x_j - x_i; % \xi
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = (norm(eta+xi) - norma)/norma; % Calculate stretch
    ee = (eta + xi)/norm(eta+xi); % Versor
    % Updating maximum stretch
    if exist('history','var') ~=0
        if S > history
            history = S;
        end
        S_max = history;
    else
        S_max = S;
        history = S;
    end
    % Defining non-existing parameters
    if ~exist('noFail','var')
        noFail = 0;
    end
    % Evaluating the force interaction
    mu = damageFactor(S_max,x_i,x_j,damage,noFail,model); % If noFail is true then we will always have mu as one

    f = c(1)*influenceFunction(norma,par_omega)*norma*fscalar(S*mu,damage)*ee; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi|  
end

function ff = fscalar(x,damage)
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
    S0 = [-0.98 0.95*Sc]; % S0- and S0+
    S1 = [-0.99 1.05*Sc]; % S1- and S1+
    % Evaluation
    if x > S1(1) && x < S0(1) % (S1-,S0-)
      ff = S0(1)*(x-S1(1))/(S0(1) - S1(1));  
    elseif x >= S0(1) && x <= S0(2) % [S0-,S0+]
      ff = x;
    elseif x > S0(2) && x < S1(2) % (S0+,S1) 
      ff = S0(2)*(S1(2)-x)/(S1(2) - S0(2));
    else
      ff = 0;
    end
else
    ff = x;
end
end
