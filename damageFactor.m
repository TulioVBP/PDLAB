function mu = damageFactor(x,x_i,x_j,noFail)
%% Input
% x: stretch, js integral, jtheta integral ...
% x_i: position of node i
% x_j: position of node j
% noFail: true if one of the nodes is in a non fail zone
% dt: timestep (not always needed)
%% Output
% mu: damage factor
%% CODE
    global damageOn crackIn model
    crackSegments = size(crackIn,1); % At least 2
    check = [];
    if damageOn
        % With damage
        for ii = 1:crackSegments-1
            [~,check(ii)] = checkBondCrack(x_i,x_j,crackIn(ii,:),crackIn(ii+1,:));
        end
        if ~sum(check)
            % Bond doesn't interceptate the initial crack (notch)
            if exist('dt','var')~=1
                dt = 0;
            end
            mu = mu_model(x);
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

function mu = mu_model(x)
% INPUT
% x - Stretch
global model S0 S1

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
        HH = @(x) (x >= 0); % Heaviside function
        xc = 0.15;
        h = @(x) 1 + HH(x).*(exp(1-1./(1-(x/xc).^2.01)) - 1) + HH(x-xc).*(-exp(1-1./(1-(x/xc).^2.01)));
        % Evaluating Ht or Hd
        mu(1) = h(x(1)); % Ht
        mu(2) = h(x(2)); % Hd-x
        mu(3) = h(x(3)); % Hd-y        
    otherwise
        mu = 1;
end

end