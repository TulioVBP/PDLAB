function mu = damageFactor(x,x_i,x_j,damage,noFail,model)
%% Input
% x: stretch, js integral, jtheta integral ...
% x_i: position of node i
% x_j: position of node j
% damage: all needed data
% noFail: true if one of the nodes is in a non fail zone
% dt: timestep (not always needed)
%% Output
% mu: damage factor
%% CODE
    crackSegments = size(damage.crackIn,1); % At least 2
    check = [];
    if damage.damageOn
        % With damage
        for ii = 1:crackSegments-1
            [~,check(ii)] = checkBondCrack(x_i,x_j,damage.crackIn(ii,:),damage.crackIn(ii+1,:));
        end
        if ~sum(check)
            % Bond doesn't interceptate the initial crack (notch)
            mu = mu_model(x,damage,model);
        else
            % Bond interceptate the initial crack (notch)
           switch model.name
           case "Lipton Free Damage"
                mu = [0 0 0];
           otherwise              
                mu = 0;
            end
        end
    else
        % No damage
        switch model.name
            case "Lipton Free Damage"
                mu = [1 1 1];
            otherwise              
                mu = 1;
        end
    end
    if exist('noFail','var')~=0
        if noFail == true
            % No damage
            switch model.name
                case "Lipton Free Damage"
                    mu = [1 1 1];
                otherwise              
                    mu = 1;
            end
        end
    end
end

function mu = mu_model(x,damage,model)
% INPUT
% x - Stretch
S0 = damage.S0;
S1 = damage.S1;

switch model.name
    case "PMB"
        if x < S0(2) && x > -1
            mu = 1;
        elseif x > S0(2) && x < S1(2)
            mu = (S1(2) - x)/(S1(2) - S0(2));
        else
            mu = 0;
        end
           
    case "Lipton Free Damage"
        %% Ht
        % Evaluating h
        xc = 0.15;
        % Evaluating Ht or Hd
        mu(1) = h_function(x(1),xc); % Ht
        mu(2) = h_function(x(2),xc); % Hd-x
        mu(3) = h_function(x(3),xc); % Hd-y        
    otherwise
        mu = 1;
end

end

function h = h_function(x,xc)
    h = 1 + heaviside_v2(x).*(exp(1-1./(1-(x/xc).^2.01)) - 1) + heaviside_v2(x-xc).*(-exp(1-1./(1-(x/xc).^2.01)));
end

function HH = heaviside_v2(x)
    HH = x >= 0; % Heaviside function
end