classdef modelPMB
    properties
        linearity = false;
        stiffnessAnal = false;
        dilatation = false;
        number = 2;
        damage;
        c;
        history;
    end
    % Methods
    methods
        function obj = modelPMB(E,par_omega,damage,G0)
            %% Pre Initialization %%
             % Any code not using output argument (obj)
             if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                E = 1;
                par_omega = [1 1 1];
                damage.DD = false;
                G0 = 1;
             end
         
             %% Object Initialization %%
             % Call superclass constructor before accessing object
             % You cannot conditionalize this statement
             %obj = obj@BaseClass1(args{:});
             %obj = modelLPST;

             %% Post Initialization %%
             % Any code, including access to object
            horizon = par_omega(1);
            mm = weightedVolume(par_omega);
            obj.c = 6*E/mm;
            if par_omega(2) == 3 && par_omega(3) == 1
                obj.damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
            elseif par_omega(2) == 1 && par_omega(3) == 0
                l = 3;
                obj.damage.Sc = sqrt((1+1/3)*pi*G0*l/(8*E*horizon*0.66467));
            elseif par_omega(2) == 1 && par_omega(3) == 1
                l = 3;
                obj.damage.Sc = sqrt(G0*(1/3+1)*pi^(3/2)/8/E*(l/horizon));
            else
                warning('Critical bond not defined.')
            end
            % Damage dependent Sc
            if damage.DD
                obj.damage.alfa = 0.2; obj.damage.beta = 0.2; obj.damage.gamma = 1.4;
            else
                obj.damage.alfa = 0; obj.damage.beta = 0; obj.damage.gamma = 1; % No dependency
            end
             
        end
        
        function [f,history,mu] = T(obj,x,u,ii,jj,dof_vec,par_omega,separatorDamage,damage,history,noFail)
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
            u = u'; % Nx1 -> 1xN
            u_i = u(dofi); u_j = u(dofj);
            xi = x_j - x_i; % \xi
            eta = u_j - u_i; % \eta
            norma = vecnorm(xi')'; 
            S = (vecnorm(eta'+xi')' - norma)./norma; % Calculate stretch
            ee = (eta + xi)./vecnorm(eta'+xi')'; % Versor
            % Updating maximum stretch
            history = obj.updateHistory(S,history);
            mu = obj.damageFactor(history',ii,1:length(jj),damage,noFail);
            %% Defining fscalar
            if damage.damageOn
                % Damage dependent crack
                alfa = obj.damage.alfa; beta = obj.damage.beta; gamma = obj.damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = obj.damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                else
                    Sc = obj.damage.Sc;
                end

                % NEW FORMULATION
                ff = (S < Sc).*S;
                % Adding noFail condition
                ff(noFail) = S(noFail);
            else
                ff = S;
            end
            
            %% Evaluating the force interaction
            f = obj.c*influenceFunction(norma,par_omega).*norma.*ff.*ee.*mu; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi| 

        end
        
        function W = strainEnergyDensity(obj,x,u,family,partialAreas,surfaceCorrection,ii,idb,par_omega,damage,historyS)
            familySet = family(family~=0);
            dofi = [idb(2*ii-1) idb(2*ii)];
            neigh_ind = 1:length(familySet);
            jj = familySet;
            xi = x(jj,:) - x(ii,:); % Nx2
            dofj = [idb(2*jj-1) idb(2*jj)]; % NX2
            u = u';
            eta = u(dofj) - u(dofi); 
            norma = vecnorm(xi,2,2);
            s = (vecnorm(xi+eta,2,2) - norma)./norma;
            noFail = damage.noFail(ii) | damage.noFail(jj);
            mu = obj.damageFactor(historyS',ii,1:length(jj),damage,noFail);
            p = mu.*s.^2/2;
            w = 1/2*obj.c*influenceFunction(norma,par_omega).*norma.^2.*p.*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)');
        end
        
        function history = updateHistory(obj,S,history)
            if nargin < 2
                error('No history was passed as argument.')
            end
            % Updating maximum stretch
            S_max = history';
            S_max(S>S_max) = S(S>S_max);
            history = S_max';
        end
        
        function mu = damageFactor(obj,x,ii,neighIndex,damage,noFail)
            %% Input
            % x: stretch, js integral, jtheta integral ...
            % ii: node i
            % neighIndex: nodes j
            % damage: all needed data
            % noFail: true if one of the nodes is in a non fail zone
            %% Output
            % mu: damage factor
            %% CODE
            if nargin <4
                error('Damage struct was not provided')
            elseif nargin < 5
                noFail = [];
            end
            
            % {Preallocate the damage factor}
            mu = ones(length(neighIndex),1);
            if damage.damageOn
                % Damage dependent crack
                alfa = obj.damage.alfa; beta = obj.damage.beta; gamma = obj.damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = obj.damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                else
                    Sc = obj.damage.Sc;
                end
                
                mu = (x< Sc)*1;
            end
            % Deal with initial broken bonds
            brokenBonds = damage.brokenBonds(ii,neighIndex);
            if any(brokenBonds)
                mu(brokenBonds,:) = zeros(sum(brokenBonds),size(mu,2));
            end
            if ~isempty(noFail)
                mu(noFail,:) = ones(size(mu(noFail,:)));
            end
        end
        
        function obj = set.history(obj,values)
            if isempty(obj.history)
                % Initialization
                x = values{1};
                maxNeigh = values{2};
                obj.history.S = zeros(length(x),maxNeigh); % S_max
            else
                % Update
                obj.history.S = values;
            end
        end
   end
end

function p = antiderivativeDTT(obj,x,damage,noFail,ii)
    % Modified PMB model
    if damage.damageOn
        % Damage dependent crack
        alfa = obj.damage.alfa; beta = obj.damage.beta; gamma = obj.damage.gamma;
        if damage.phi(ii) > alfa
            Sc = obj.damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
        else
            Sc = obj.damage.Sc;
        end
        S0 = [-0.98 0.95*Sc]; % S0- and S0+
        S1 = [-0.99 1.05*Sc]; % S1- and S1+
        % Evaluate integration constants
        A = [1 0 0 0 0; 0 1 0 0 0; 1 0 0 -1 0; 0 0 1 0 0; 0 0 -1 0 1];
        b = [S0(1)^2/2 - S0(1)/(S0(1) - S1(1))*(S0(1)^2/2 - S1(1)*S0(1));
            0;
            S0(1)/(S0(1) - S1(1))*(S1(1)^2/2);
            S0(2)^2/2 - S0(2)/(S1(2) - S0(2))*(-S0(2)^2/2 + S1(2)*S0(2));
            S0(2)/(S1(2) - S0(2))*S1(2)^2/2];
        C = A\b;
        p = (x<=S1(1)).* C(4) + (x<=S0(1)).*(x>S1(1)).*(S0(1)/(S0(1) - S1(1)).*(x.^2/2 - S1(1)*x) + C(1)) ...
            + (x<=S0(2)).*(x>S0(1)).*(x.^2/2 + C(2)) + (x<=S1(2)).*(x>S0(2)).*(S0(2)/(S1(2) - S0(2)).*(S1(2)*x - x.^2/2) + C(3)) ...
            + (x>S1(2)).*C(5);
        % {Correcting the noFail}
        p(noFail) = x(noFail).^2/2; 
    else
        p = x.^2/2;
    end
end