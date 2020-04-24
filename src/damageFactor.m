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
    if model.number == 3 && model.number == 6 % Lipton
        mu = zeros(length(neighIndex),3);
    else
        mu = zeros(length(neighIndex),1);
    end
    
    if damage.damageOn
        brokenBonds = damage.brokenBonds(ii,neighIndex);
        % MU MODEL
         switch model.number
            case 1 %"PMB DTT"
                % Damage dependent crack
                alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                else
                    Sc = damage.Sc;
                end
                S0 = [-0.98 0.95*Sc]; % S0- and S0+
                S1 = [-0.99 1.05*Sc]; % S1- and S1+
                mu = (x<= S0(2)).*(x>=-1).*1 + (x > S0(2)).*(x<S1(2)).*(S1(2) - x)/(S1(2) - S0(2));
            case 3 %"Lipton Free Damage"
               %% Ht
               % Evaluating h
               %xc = 0.15;
               %xc = 0.15*10^-6;
               xc = (0.05)^2/(1+1.05^2) * 0.02e-6; % js(Sc)*dt = 2.3781e-11
               mu = (x<xc).*(exp(1-1./(1-(x/xc).^2.01)));
               mu(isnan(mu)) = zeros(sum(sum(isnan(mu))),1);
            case 6 %"LSJ-T"
               %% Ht
               % Evaluating h
               %xc = 0.15*10^-6;
               xc = (0.05)^2/(1+1.05^2) * 0.02e-6; % js(Sc)*dt = 2.3781e-11
               mu = (x<xc).*(exp(1-1./(1-(x/xc).^2.01)));
               mu(isnan(mu)) = zeros(sum(sum(isnan(mu))),1);
            case 4 % "LPS 2D"
                 % Damage dependent crack
                 alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                else
                    Sc = damage.Sc;
                end
                S0 = [-0.98 0.95*Sc]; % S0- and S0+
                S1 = [-0.99 1.05*Sc]; % S1- and S1+
                mu = (x<= S0(2)).*(x>=-1).*1 + (x > S0(2)).*(x<S1(2)).*(S1(2) - x)/(S1(2) - S0(2));
            case 5 % "PMB"
                % Damage dependent crack
                alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                    if length(Sc)>1
                        disp('error');
                    end
                else
                    Sc = damage.Sc;
                end
                mu = (x<Sc) * 1;
            otherwise
                mu = 1;
        end
        if any(brokenBonds)  % Bond interceptate the initial crack (notch)
            switch model.number
               case 3 %"Lipton Free Damage"
                    mu(brokenBonds,1) = zeros(sum(brokenBonds),1);
               case 6 %"LSJ-T"
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

switch model.number
    case 1 %"PMB-DTT"
        % Damage dependent crack
        alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
        if damage.phi > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
        else
            Sc = damage.Sc;
        end
        S0 = [-0.98 0.95*Sc]; % S0- and S0+
        S1 = [-0.99 1.05*Sc]; % S1- and S1+
        mu = (x<= S0(2)).*(x>=-1).*1 + (x > S0(2)).*(x<S1(2)).*(S1(2) - x)/(S1(2) - S0(2));
    case 3 %"Lipton Free Damage"
        %% Ht
        % Evaluating h
        xc = 0.15*10^-6;
        % Evaluating Ht or Hd
       mu = (x<xc).*(exp(1-1./(1-(x/xc).^2.01)));
       mu(isnan(mu)) = zeros(sum(sum(isnan(mu))),1);
    case 5 % "PMB"
        % Damage dependent crack
        alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
        if damage.phi > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi-alfa)/(1-damage.phi));
        else
            Sc = damage.Sc;
        end
        mu = (x<Sc) * 1;
    case 7 % "PMB Concrete"
        Sc = [damage.Sc damage.St];
        mu = (x<Sc(2) & x > Sc(1))*1;
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
