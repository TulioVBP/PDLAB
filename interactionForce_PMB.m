function [f,S_max] = interactionForce_PMB(x_i,x_j,u_i,u_j,S_max_ant,notch,noFail)
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
    xi = x_j - x_i; % \xi
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = (norm(eta+xi) - norm(xi))/norma; % Calculate stretch
    ee = (eta + xi)/norm(eta+xi); % Versor
    % Updating maximum stretch
    if S > S_max_ant
        S_max = S;
    else
        S_max = S_max_ant;
    end
    % Evaluating the force interaction
    f = c1*influenceFunction(norma,horizon,omega)*norma*fscalar(S*damageFactor(S_max,notch,x_i,x_j)*noFail)*ee; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi|   
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

%function mu = damageFactor(x)
%     global S0 S1
%     if x < S0(2) && x > -1
%         mu = 1;
%     elseif x > S0(2) && x < S1(2)
%         mu = (S1(2) - x)/(S1(2) - S0(2));
%     else
%         mu = 0;
%     end
% end