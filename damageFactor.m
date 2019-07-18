function mu = damageFactor(x,ii,neighIndex,damage,noFail,model)
%% Input
% x: stretch, js integral, jtheta integral ...
% x_i: position of node i
% x_j: position of node j
% damage: all needed data
% noFail: true if one of the nodes is in a non fail zone
% model: timestep (not always needed)
% brokenBonds: 
%% Output
% mu: damage factor
%% CODE
    crackSegments = size(damage.crackIn,1); % At least 2
    check = zeros(length(neighIndex),crackSegments-1);
    % {Preallocate the damage factor}
    if model.dilatation
        mu = zeros(length(neighIndex),3);
    else
        mu = zeros(length(neighIndex),1);
    end
    
    if damage.damageOn
        brokenBonds = damage.brokenBonds(ii,neighIndex);
        mu = mu_model(x,damage,model);
        if any(brokenBonds)  % Bond interceptate the initial crack (notch)
            switch model.name
                case "Lipton Free Damage"
                    mu(brokenBonds,1) = zeros(sum(brokenBonds),1);
               otherwise              
                    mu(brokenBonds,:) = zeros(sum(brokenBonds),size(mu,2));
            end
        end
    else
        % No damage
        mu = ones(size(mu));
    end
    if ~isempty(noFail)%exist('noFail','var')~=0
        mu(noFail,:) = ones(size(mu(noFail,:)));
    end
end

function mu = mu_model(x,damage,model)
% INPUT
% x - Stretch

switch model.name
    case "PMB"
        % Damage dependent crack
        alfa = 0.2; beta = 0.2; gamma = 1.4;
        if damage.phi > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
        else
            Sc = damage.Sc;
        end
        S0 = [-0.98 0.95*Sc]; % S0- and S0+
        S1 = [-0.99 1.05*Sc]; % S1- and S1+
%         if x < S0(2) && x > -1
%             mu = 1;
%         elseif x > S0(2) && x < S1(2)
%             mu = (S1(2) - x)/(S1(2) - S0(2));
%         else
%             mu = 0;
%         end
        mu = (x<= S0(2)).*(x>=-1).*1 + (x > S0(2)).*(x<S1(2)).*(S1(2) - x)/(S1(2) - S0(2));
    case "Lipton Free Damage"
        %% Ht
        % Evaluating h
        xc = 0.15*10^-6;
        % Evaluating Ht or Hd
        mu(:,1) = h_function(x(:,1),xc); % Ht
        mu(:,2) = h_function(x(:,2),xc); % Hd-x
        mu(:,3) = h_function(x(:,3),xc); % Hd-y        
    otherwise
        mu = 1;
end

end

function h = h_function(x,xc)
%     if x < xc
        h = (x<xc).*(exp(1-1./(1-(x/xc).^2.01)));
        if any(isnan(h))
           h(isnan(h)) = zeros(sum(isnan(h)),1);
        end
%     else
%         h = 0;
%     end
end
