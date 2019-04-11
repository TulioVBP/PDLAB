function [f,S_max] = interactionForce_LLPSBB(x_i,x_j,u_i,u_j,S_max_ant,notch)
%% Function to evaluate the linearized LPS bond-based force between two interacting nodes
%% INPUT
% x_i: position of node i
% x_j: position of node j
% u_i: displacement of node i
% u_j: displacement of node j
% S_max_ant: maximum stretch for each given bond
% notch: coordinates of the initial notch
%% OUTPUT
% f: vector state force between j and i nodes
% S_max: maximum stretch for each bond
%% CODE
    global c1 horizon omega
    xi = x_j - x_i; % \xi
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = (norm(eta+xi) - norm(xi))/norma; % Calculate stretch
    ee = xi/norma; % Versor
    % Updating maximum stretch
    if S > S_max_ant
        S_max = S;
    else
        S_max = S_max_ant;
    end
    % Evaluating the force interaction
    f = c1*influenceFunction(norma,horizon,omega)*fscalar(eta,ee,S)*ee;    
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
