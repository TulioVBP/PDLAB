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
    u = u';
    u_i = u(dofi); u_j = u(dofj);
    eta = u_j - u_i; % \eta
    norma = vecnorm(xi')';
    %xi_rep = repmat(xi,size(x_j,1),1);
    S = dot(eta',xi')'./norma.^2; % Calculate stretch - linear
    ee = xi./norma; % Versor
    % Updating maximum stretch
    if nargin > 10  && damage.damageOn% Damage considered
        S_max = history';
        history(S>S_max) = S(S>S_max);
        S_max = history;
        % Evaluating the damage factor
        mu = damageFactor(S_max,ii,1:length(jj),damage,noFail,model); % If noFail is true then we will always have mu as one
    else % No damage considered
        history = zeros(length(S),1);
        mu = ones(length(S),1);
        noFail = [];
    end
%     if ~exist('damage','var')
%         damage.damageOn = false;
%         damage.crackIn = [];
%     end
    f = c(1)*influenceFunction(norma,par_omega).*fscalar(eta,ee,S,damage).*ee;
end

function ff = fscalar(eta,versor,x,damage)
% eta: uj - ui
% versor: direction of the force
% x: stretch
if damage.damageOn
        %versor_rep = repmat(versor,size(eta,1),1);
        ff = (x < damage.S0(2)).*dot(eta',versor')';
else
    ff = dot(eta',versor')';
end
end
