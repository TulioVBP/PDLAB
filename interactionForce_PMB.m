function [f,history] = interactionForce_PMB(x,u,ii,dof_vec,familyMat,neighIndex,dt,history,noFail)
%% INPUT
% x_i: position of node i
% x_j: position of node j
% u_i: displacement of node i
% u_j: displacement of node j
% S_max_ant: maximum stretch for each given bond
% notch: coordinates of the initial notch
% noFail: true if the damage is off for this specific bond
%% OUTPUT
% f: vector state force between j and i nodes
% S_max: maximum stretch for each bond
%% CODE
    global c1 horizon omega
    jj = familyMat(ii,neighIndex);
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
            S_max = S;
        else
            S_max = history;
        end
    else
        S_max = S;
        history = S;
    end
    % Defining non-existing parameters
    if ~exist('noFail','var')
        noFail = 0;
    end
    % Evaluating the force interaction
    mu = damageFactor(S_max,x_i,x_j,noFail); % If noFail is true then we will always have mu as one
    f = c1*influenceFunction(norma,horizon,omega)*norma*fscalar(S*mu)*ee; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi|  
end

function ff = fscalar(x)
global S0 S1 damageOn
if damageOn
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
