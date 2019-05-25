function [f,history,mu] = interactionForce_LLPSBB(x,u,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
%% Function to evaluate the linearized LPS bond-based force between two interacting nodes
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
    xi = x_j - x_i; % \xi
    u_i = u(dofi)'; u_j = u(dofj)';
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = dot(eta,xi)/norma^2; % Calculate stretch - linear
    ee = xi/norma; % Versor
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
    % Evaluating the force interaction
    if exist('noFail','var') == 0
        noFail = 1;
    end
    if ~exist('damage','var')
        damage.damageOn = false;
        damage.crackIn = [];
    end
    f = c(1)*influenceFunction(norma,par_omega)*fscalar(eta,ee,S,damage)*ee;
    history = 0; % For this specific model
    mu = 1; % No damage in this model
end

function ff = fscalar(eta,versor,x,damage)
% eta: uj - ui
% versor: direction of the force
% x: stretch
if damage.damageOn
     if x < damage.S0(2)
        ff = dot(eta,versor);
     else
        ff = 0;
     end
else
    ff = dot(eta,versor);
end
end
