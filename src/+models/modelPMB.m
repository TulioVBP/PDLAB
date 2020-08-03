classdef modelPMB
    properties
        b_linearity = false;
        b_stiffnessAnal = true;
        b_dilatation = false;
        number = 2;
        damage;
        c;
        history;
    end
    % Methods
    methods
        %% CONSTRUCTOR
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
        %% Force vector state
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
        
        %% Strain energy density
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
        
        %% Update history
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
        
        %% Stiffness matrix
        function A = analyticalStiffnessMatrix(obj,x,u,ndof,idb,familySet,partialAreas,surfaceCorrection,V,par_omega,damage,history,mu)
            C = obj.c(1)/2*weightedVolume(par_omega);
            u = u';
            penalty = 1e10;
            N = size(x,1);
            A = zeros(2*N,2*N); % 2N GDLs
            m = weightedVolume(par_omega)*ones(length(x),1);
           
            for ii = 1:N
                dofi = [idb(2*ii-1) idb(2*ii)];
                family = familySet(ii,familySet(ii,:)~=0);
                iII = 1:length(family);
                jj = family'; % j sum
                dofj = [idb(2*jj-1) idb(2*jj)];
                eta = u(dofj) - u(dofi);
                xi = x(jj,:) - x(ii,:);
                M = eta+xi;
                normaj = vecnorm(xi')'; 
                norma_eta = vecnorm(M')';
                omegaj = influenceFunction(normaj,par_omega);
                muj = 1;
                if nargin > 12 % damage factor provided
                    muj = mu{ii};
                elseif nargin > 10 % Damage parameters provided
                   S = (vecnorm(eta+xi,2,2) - normaj)./normaj; % Calculate stretch
                   history_S = obj.updateHistory(S,history.S(ii,iII));
                   noFail = damage.noFail(ii) | damage.noFail(jj);
                   muj = obj.damageFactor(history_S',ii,1:length(jj),damage,noFail); 
                end               
                % U
                if dofi(1) <= ndof
                    % First dof of node ii is free
                    ti1u = C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aii
                    ti2u = C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aip
                    tj1u = -C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aij
                    tj2u = -C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aijp
                    A(dofi(1),dofi(1)) = sum(ti1u);
                    A(dofi(1),dofj(:,1)) = tj1u';
                    A(dofi(1),dofi(2)) = sum(ti2u);
                    A(dofi(1),dofj(:,2)) = tj2u';
               else
                    % Constraint nodes
                    A(dofi(1),dofi(1)) = -penalty;
               end
               if dofi(2) <= ndof
                  % V
                   ti1v = C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   ti2v = C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   tj1v = -C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   tj2v = -C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
                   A(dofi(2),dofi(1)) = sum(ti1v);
                   A(dofi(2),dofj(:,1)) = tj1v';
                   A(dofi(2),dofi(2)) =sum(ti2v);
                   A(dofi(2),dofj(:,2)) = tj2v';
               else
                    % Constraint nodes
                    A(dofi(2),dofi(2)) = -penalty;
               end
            end
            A = -A;
        end
   end
end

