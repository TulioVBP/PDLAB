function [f,S_max] = interactionForce(x_i,x_j,u_i,u_j,S_max_ant)
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
    f = c1*influenceFunction(norma,horizon,omega)*fscalar(S*damageFactor(S_max))*ee;    
end

function ff = fscalar(x)
global S0 S1
if x > S1(1) && x < S0(1) % (S1-,S0-)
  ff = S0(1)*(x-S1(1))/(S0(1) - S1(1));  
elseif x >= S0(1) && x <= S0(2) % [S0-,S0+]
  ff = x;
elseif x > S0(2) && x < S1(2) % (S0+,S1) 
  ff = S0(2)*(S1(2)-x)/(S1(2) - S0(2));
else
  ff = 0;
end
end

function mu = damageFactor(x)
    global S0 S1
    if x < S0(2) && x > -1
        mu = 1;
    elseif x > S0(2) && x < S1(2)
        mu = (S1(2) - x)/(S1(2) - S0(2));
    else
        mu = 0;
    end
end