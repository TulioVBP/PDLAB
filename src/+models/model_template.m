classdef model_template
    properties
        %% Mandatory properties
        b_linearity = false; % True if the model is linear with the displacement field
        b_stiffnessAnal = true; % True if there is a analyticalStiffnessMatrix implemented
        b_dilatation = false; % True if the model has a dilatation term
        history; % Property that accumulates the history of a specific variable (e.g., maximum bond stretch for PMB)
        %% Recommended properties
        damage; % Struct property helpful for damage model
        c; % Material constants
    end
    %% Methods
    methods
        %% CONSTRUCTOR METHOD
        function obj = model_template(E,par_omega,damage,G0)
            % Pre Initialization %%
             % Any code not using output argument (obj)
             if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                E = 1; % Young modulus
                par_omega = [1 1 1]; % [horizon omega_option omega_gamma]
                G0 = 1; % Energy release rate
             end

            % Post Initialization %%
            % Any code, including access to object
            horizon = par_omega(1);
            mm = weightedVolume(par_omega);
            % Random initialization
            obj.c = E;
            obj.damage.Sc = G0/4;
             
        end
        %% Force vector state
        function [f,history,mu] = T(obj,x,u,ii,jj,dof_vec,par_omega,separatorDamage,damage,history,noFail)
            % ----- INPUT
            % x - node position matrix
            % u - degree of freedom displacement vector
            % ii - index of the i-th node
            % jj - index of the j-th node inside i's family
            % dof_vec - matrix with the degrees of freedom corresponding for each node
            % separatorDamage - doesn't do anything but is useful to separate the
            %                   normal variables to the ones needed for damage simulation
            % damage - struct with useful info for damage
            % history - maximum stretch for the given bond
            % nofail - boolean variable that takes true if 
            % ----- OUTPUT
            % f: force interaction between j and i nodes. In state-based: f
            % =  T(x)<x-x'> - T(x')<x'-x>
            % history: maximum stretch for each bond
            % mu: damage factor used in the damage index variable
            %% CODE
            x_i = x(ii,:); % Size n x 1
            x_j = x(jj,:); % Size: n x 2
            dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
            u = u'; % Nx1 -> 1xN
            u_i = u(dofi); % Size: 1 x 2
            u_j = u(dofj); % Size: n x 2
            xi = x_j - x_i; % Size: n x 2
            eta = u_j - u_i; % Size: n x 2
            norma = vecnorm(xi')'; 
            S = (vecnorm(eta'+xi')' - norma)./norma; % Calculate stretch
            ee = (eta + xi)./vecnorm(eta'+xi')'; % Versor
            
            % Updating maximum stretch
            history = obj.updateHistory(S,history);
            mu = obj.damageFactor(history',ii,1:length(jj),damage,noFail);
            
            %% Evaluating the force interaction
            f = obj.c*influenceFunction(norma,par_omega).*norma.*ff.*ee.*mu;

        end
        
        %% Strain energy density
        function W = strainEnergyDensity(obj,x,u,family,partialAreas,surfaceCorrection,ii,idb,par_omega,damage,historyS)
            % ----- INPUT
            % x - node position matrix
            % u - degree of freedom displacement vector
            % family - set of neighbors, so it can be [j1, j2, ..., jn, 0, ... , 0]
            % partialAreas - partial areas of the corresponding neighbors
            % surfaceCorrection - surface effects correction factor
            % ii - node at where we are evaluating the strain energy
            % idb - every line corresponds to the global degree of freedom
            %       of the node ii (i.e., lines 2 ii -1 and 2 ii
            %       corresponds to the dof of node ii), and the assigned
            %       value is the corresponding rearranged DOF
            % par_omega - [horizon omega_option gamma]
            % damage - damage struct
            % historyS - history variable for S
            % ------ OUTPUT
            % W - strain energy density at x_ii
            
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
            % Updating maximum stretch - Example, passing the maximum
            % stretch
            S_max = history';
            S_max(S>S_max) = S(S>S_max);
            history = S_max';
        end
        
        %% Damage factor
        function mu = damageFactor(obj,x,ii,neighIndex,damage,noFail)
            % ----- Input
            % x: stretch
            % ii: node i
            % neighIndex: nodes j
            % damage: all needed data
            % noFail: true if one of the nodes is in a non fail zone
            % ----- Output
            % mu: damage factor
            %% CODE
            if nargin <4
                error('Damage struct was not provided')
            elseif nargin < 5
                noFail = [];
            end
            
            % {Preallocate the damage factor}
            mu = (x< Sc)*1;
            
            % Deal with initial broken bonds
            brokenBonds = damage.brokenBonds(ii,neighIndex);
            if any(brokenBonds)
                mu(brokenBonds,:) = zeros(sum(brokenBonds),size(mu,2));
            end
            if ~isempty(noFail)
                mu(noFail,:) = ones(size(mu(noFail,:)));
            end
        end
        
        %% Initialization of the history variable
        function obj = set.history(obj,values)
            if isempty(obj.history)
                % Initialization - For the first call in experiments
                % template
                x = values{1};
                maxNeigh = values{2};
                obj.history.S = zeros(length(x),maxNeigh); % S_max
            else
                % Update - In case one wants to update the objects history
                obj.history.S = values;
            end
        end
        
        %% Stiffness matrix
        function A = analyticalStiffnessMatrix(obj,x,u,ndof,idb,familySet,partialAreas,surfaceCorrection,V,par_omega,damage,history)
            % ----- INPUT
            % x - node position matrix
            % u - degree of freedom displacement vector
            % ndof - number of degrees of freedom
            % familySet - set of neighbors, so it can be [j1, j2, ..., jn, 0, ... , 0]
            % partialAreas - partial areas of the corresponding neighbors
            % surfaceCorrection - surface effects correction factor
            % V - vector with cell's area (e.g. h^2 for homogeneous grid)
            % par_omega - [horizon omega_option gamma]
            % damage - damage struct
            % history - history variable (it can be the struct)
            
            % ------ OUTPUT
            % W - strain energy density at x_ii
            
            C = obj.c(1)/2*weightedVolume(par_omega);
            u = u';
            penalty = 1e10;
            N = size(x,1);
            A = zeros(2*N,2*N); % 2N GDLs
            m = weightedVolume(par_omega)*ones(length(x),1);
            
            % EXAMPLE FOR PMB
%             for ii = 1:N
%                 dofi = [idb(2*ii-1) idb(2*ii)];
%                 family = familySet(ii,familySet(ii,:)~=0);
%                 iII = 1:length(family);
%                 jj = family'; % j sum
%                 dofj = [idb(2*jj-1) idb(2*jj)];
%                 eta = u(dofj) - u(dofi);
%                 xi = x(jj,:) - x(ii,:);
%                 M = eta+xi;
%                 normaj = vecnorm(xi')'; 
%                 norma_eta = vecnorm(M')';
%                 omegaj = influenceFunction(normaj,par_omega);
%                 muj = 1;
%                 if nargin > 10 % Damage parameters provided
%                    S = (vecnorm(eta+xi,2,2) - normaj)./normaj; % Calculate stretch
%                    history_S = obj.updateHistory(S,history.S(ii,iII));
%                    noFail = damage.noFail(ii) | damage.noFail(jj);
%                    muj = obj.damageFactor(history_S',ii,1:length(jj),damage,noFail); 
%                 end               
%                 % U
%                 if dofi(1) <= ndof
%                     % First dof of node ii is free
%                     ti1u = C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aii
%                     ti2u = C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii); % Aip
%                     tj1u = -C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aij
%                     tj2u = -C*(1./m(jj) + 1/m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,1).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);% Aijp
%                     A(dofi(1),dofi(1)) = sum(ti1u);
%                     A(dofi(1),dofj(:,1)) = tj1u';
%                     A(dofi(1),dofi(2)) = sum(ti2u);
%                     A(dofi(1),dofj(:,2)) = tj2u';
%                else
%                     % Constraint nodes
%                     A(dofi(1),dofi(1)) = -penalty;
%                end
%                if dofi(2) <= ndof
%                   % V
%                    ti1v = C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
%                    ti2v = C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
%                    tj1v = -C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,1)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
%                    tj2v = -C*(1./m(jj) + 1./m(ii)).*(omegaj./(norma_eta).^2.*M(:,2)).*M(:,2).*partialAreas(ii,iII)'.*surfaceCorrection(ii,iII)'.*muj.*V(ii);
%                    A(dofi(2),dofi(1)) = sum(ti1v);
%                    A(dofi(2),dofj(:,1)) = tj1v';
%                    A(dofi(2),dofi(2)) =sum(ti2v);
%                    A(dofi(2),dofj(:,2)) = tj2v';
%                else
%                     % Constraint nodes
%                     A(dofi(2),dofi(2)) = -penalty;
%                end
%             end
%             A = -A;
        end
   end
end

