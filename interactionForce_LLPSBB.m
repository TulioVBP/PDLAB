function [f,history] = interactionForce_LLPSBB(x,u,ii,dof_vec,familyMat,partialAreas,neighIndex,dt,history,noFail)
%% Function to evaluate the linearized LPS bond-based force between two interacting nodes
%% INPUT
% x_i: position of node i
% x_j: position of node j
% u_i: displacement of node i
% u_j: displacement of node j
% S_max_ant: maximum stretch for each given bond
% notch: coordinates of the initial notch
% noFail: true if the damage is off for this specific bond
%% OUTPUT
% f: vector internal force acting on node i due to the j-th node on its
%    neighbourhood
% T: vector state force acting on i by j 
% S_max: maximum stretch for each bond
%% CODE
    global c1 horizon omega
    jj = familyMat(ii,neighIndex);
    x_i = x(ii,:); x_j = x(jj,:);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x_j - x_i; % \xi
    u_i = u(dofi)'; u_j = u(dofj)';
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = dot(eta,xi)/norma^2; % Calculate stretch - linear
    ee = xi/norma; % Versor
    % Updating maximum stretch
    if exist('S_max_ant','var') ~=0
        if S > S_max_ant
            S_max = S;
        else
            S_max = S_max_ant;
        end
    else
        S_max = S;
    end
    % Evaluating the force interaction
    if exist('noFail','var') == 0
        noFail = 1;
    end
    f = c1*influenceFunction(norma,horizon,omega)*fscalar(eta,ee,S)*noFail*ee;
    history = 0; % For this specific model
end

function ff = fscalar(eta,versor,x)
% eta: uj - ui
% versor: direction of the force
% x: stretch
global S0 S1 damageOn
if damageOn
     if x < S0(2)
        ff = dot(eta,versor);
     else
        ff = 0;
     end
else
    ff = dot(eta,versor);
end
end
