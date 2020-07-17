classdef modelDTT
    properties
        linearity = true;
        stiffnessAnal = true;
        dilatation = false;
        number = 1;
        damage;
        c;
    end
    % Methods
    methods
        function obj = modelDTT(E,parOmega)
            %% Pre Initialization %%
             % Any code not using output argument (obj)
             if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                args{1} = 1;
                args{2} = [1 1 1];
             else
                % When nargin ~= 0, assign to cell array,
                % which is passed to supclass constructor
                args{1} = E;
                args{2} = parOmega;
             end
         
             %% Object Initialization %%
             % Call superclass constructor before accessing object
             % You cannot conditionalize this statement
             %obj = obj@BaseClass1(args{:});
             %obj = modelLPST;

             %% Post Initialization %%
             % Any code, including access to object
             mm = weightedVolume(args{2});
             obj.c = 6*args{1}/mm;
        end
        
        function [f,history,mu] = T(x,u,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
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
            %if exist('history','var') ~=0
            if nargin > 10  && damage.damageOn% Damage considered
                S_max = history';
                history(S>S_max) = S(S>S_max);
                S_max = history';
                % Evaluating the damage factor
                mu = damageFactor(S_max,ii,1:length(jj),damage,noFail,model); % If noFail is true then we will always have mu as one
            else % No damage considered
                history = S;
                mu = ones(length(S),1);
                noFail = [];
            end
            %% Defining fscalar
            x = S.*mu;
            if damage.damageOn
                % Damage dependent crack
                alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
                if damage.phi(ii) > alfa
                    Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
                else
                    Sc = damage.Sc;
                end
                S0 = [-0.98 0.95*Sc]; % S0- and S0+
                S1 = [-0.99 1.05*Sc]; % S1- and S1+

                % NEW FORMULATION
                ff = (x>S1(1)).*(x<S0(1)).*(S0(1)*(x-S1(1))./(S0(1) - S1(1))) + (x>=S0(1)).*(x<=S0(2)).*x ...
                    + (x>S0(2)).*(x<S1(2)).*(S0(2)*(S1(2)-x)./(S1(2)-S0(2)));
                % Adding noFail condition
                ff(noFail) = x(noFail);
            else
                ff = x;
            end
            
            %% Evaluating the force interaction
            f = c(1)*influenceFunction(norma,par_omega).*norma.*ff.*ee; % Influence function times norma because the omega_d used is related to the original influence function by omega_d = omega*|\xi| 

        end
        
        function W = strainEnergyDensity(x,u,theta,family,partialAreas,surfaceCorrection,ii,idb,par_omega,c,model,damage,historyS)
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
            
            if  damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                mu = modelDTT.damageFactor(historyS(neigh_ind)',ii,neigh_ind,damage,noFail); % NoFail not required
            else
                noFail = [];
                mu = ones(length(jj),1);
            end
            p = antiderivativeDTT(s,damage,noFail,ii);
            w = 1/2*c(1)*influenceFunction(norma,par_omega).*norma.^2.*p.*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)');
        end
        
        function historyS = updateHistory(S,historyS)
            if nargin < 2
                error('No history was passed as argument.')
            end
            % Updating maximum stretch
            S_max = historyS;
            historyS(S>S_max) = S(S>S_max);
        end
        
        function mu = damageFactor(x,ii,neighIndex,damage,noFail)
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
            mu = zeros(length(neighIndex),1);
            if damage.damageOn
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
            
   end
end

function p = antiderivativeDTT(x,damage,noFail,ii)
    % Modified PMB model
    if damage.damageOn
        % Damage dependent crack
        alfa = damage.alfa; beta = damage.beta; gamma = damage.gamma;
        if damage.phi(ii) > alfa
            Sc = damage.Sc*min(gamma,1+beta*(damage.phi(ii)-alfa)/(1-damage.phi(ii)));
        else
            Sc = damage.Sc;
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