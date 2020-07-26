classdef modelLPST
     properties
        linearity = false;
        b_stiffnessAnal = true;
        b_dilatation = true;
        number = 1;
        damage;
        c;
        history;
        theta;
     end
    methods
        function obj = modelLPST(E,nu,par_omega,damage,G0)
            %% Pre Initialization %%
             % Any code not using output argument (obj)
             if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                E = 1;
                nu = 1/4;
                par_omega = [1 1 1];
                damage.DD = false;
                G0 = 0;
             end
            mm = weightedVolume(par_omega);
            kappa = E/3/(1-2*nu); 
            mu = E/2/(1+nu);
            horizon = par_omega(1);
            obj.c(1) = kappa + mu/9*(nu+1)^2/(2*nu-1)^2;
            obj.c(2) = 8*mu/mm;
            obj.c(3) = nu;
            if par_omega(2) == 3 && par_omega(3) == 1
                obj.damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E)/horizon);
            elseif par_omega(2) == 1 && par_omega(3) == 0
                l = 4.3;
                obj.damage.Sc = sqrt((1+nu)*pi*G0*l/(8*E*horizon*0.66467));
            elseif par_omega(2) == 1 && par_omega(3) == 1
                l = 3;
                obj.damage.Sc = sqrt(G0*(nu+1)*pi^(3/2)/8/E*(l/horizon));
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
        
        function [theta,historyT] = dilatation(obj,x,u,family,partialAreas,surfaceCorrection,transvList,idb,par_omega,damage,historyS,historyT)
            m = weightedVolume(par_omega);
            if isempty(transvList) % Not a specific range of nodes was chosen
               transvList = 1:length(x);
            end
            
            theta = zeros(length(transvList),1); % Initialize dilatation vector
            transv_ind = 1; % Introduced so that we can pass as argument a smaller matrix
            u = u';
            for ii = transvList
                dofi = [idb(2*ii-1) idb(2*ii)];
                familySet = family(transv_ind,family(transv_ind,:)~=0);
                jj = familySet;
                neigh_ind = 1:length(jj);
                dofj = [idb(2*jj-1) idb(2*jj)];
                xi = x(jj,:)-x(ii,:);
                norma = vecnorm(xi,2,2);
                eta = u(dofj)-u(dofi);
                nu = obj.c(3);
                elong = vecnorm(xi+eta,2,2) - vecnorm(xi,2,2);
                S = elong./norma;
                theta_vec = 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,par_omega).*norma.*elong.*partialAreas(transv_ind,neigh_ind)'.*surfaceCorrection(transv_ind,neigh_ind)';
                % Updating maximum stretch
                history_upS = obj.updateHistory(S,historyS(ii,neigh_ind)');
                noFail = damage.noFail(ii) | damage.noFail(jj);
                mu = obj.damageFactor(history_upS,ii,1:length(jj),damage,noFail);
                
                theta(transv_ind) = sum(theta_vec.*mu);
                
                transv_ind = transv_ind + 1;
            end
        end
        
        function [f,historyS,mu] = T(obj,x,u,theta,ii,jj,dof_vec,par_omega,separatorDamage,damage,historyS,historyTheta,noFail)
            nu = obj.c(3);
            dofi = dof_vec(ii,:); 
            dofj = dof_vec(jj,:);
            xi = x(jj,:) - x(ii,:); % \xi
            u = u';
            u_i = u(dofi); u_j = u(dofj);
            eta = u_j - u_i; % \eta
            norma = vecnorm(xi,2,2); 
            elong = vecnorm(xi+eta,2,2) - norma; % Calculate elongation - linear
            S = elong./norma;
            ee = (xi+eta)./vecnorm(xi+eta,2,2); % Versor
            
             % Updating maximum stretch
            historyS = obj.updateHistory(S,historyS');
            noFail = damage.noFail(ii) | damage.noFail(jj);
            mu = obj.damageFactor(historyS,ii,1:length(jj),damage,noFail);
            
            % ---- Evaluatin the force ----
            % Dilatation term
            %theta = [theta(ii) theta(jj)];
            m = weightedVolume(par_omega);
            ff = fscalar(obj,S,damage,noFail,ii);
            T_ij = 2*(2*nu-1)/(nu-1)*((obj.c(1) + obj.c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega).*mu.*norma/m)*theta(ii) ...
                + obj.c(2)*influenceFunction(norma,par_omega).*norma.*mu.*ff;%.*(elong);
            T_ji = 2*(2*nu-1)/(nu-1)*((obj.c(1) + obj.c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega).*mu.*norma/m).*theta(jj) ...
                + obj.c(2)*influenceFunction(norma,par_omega).*norma.*mu.*ff;
            f = (T_ij + T_ji).*ee;
            
        end
        
        function historyS = updateHistory(obj,S,historyS)
            if nargin < 2
                error('No history was passed as argument.')
            end
            % Updating maximum stretch
            historyS(S>historyS) = S(S>historyS);
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
        
        function W = strainEnergyDensity(obj,x,u,theta,family,partialAreas,surfaceCorrection,ii,idb,par_omega,damage,historyS,historyT) 
            familySet = family(family~=0);
            dofi = [idb(2*ii-1) idb(2*ii)];
            m = weightedVolume(par_omega);
            horizon = par_omega(1);
            % Evaluate dilatation
            theta_i = theta(ii);
            neigh_ind = 1:length(familySet);
            jj = familySet;
            xi = x(jj,:) - x(ii,:); % Nx2
            dofj = [idb(2*jj-1) idb(2*jj)]; % NX2
            u = u';
            eta = u(dofj) - u(dofi); 
            norma = vecnorm(xi')';
            s = (vecnorm(xi'+eta')' - norma)./norma;
            
            noFail = damage.noFail(ii) | damage.noFail(jj);
            mu = obj.damageFactor(historyS(neigh_ind)',ii,neigh_ind,damage,noFail); % NoFail not required
            
            p = antiderivativeDTT(obj,s,damage,noFail,ii);
            nu = obj.c(3);
            w = obj.c(2)/2*influenceFunction(norma,par_omega).*norma.^2.*(2*p).*mu;
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)')+ (obj.c(1)/2 + obj.c(2)*m/3*(1/6 - (nu-1)/(2*(2*nu-1))))*theta_i^2;
        end
        
        function obj = set.history(obj,values)
            if iscell(values)
                % Initialization
                x = values{1};
                maxNeigh = values{2};
                obj.history.S = zeros(length(x),maxNeigh); % S_max
                obj.history.theta = zeros(length(x),1);
            else
                % Update
                obj.history.S = values;
            end
        end
        
        function A = analyticalStiffnessMatrix(obj,x,u,ndof,idb,familySet,partialAreas,surfaceCorrection,V,par_omega,damage,history,mu)
            u = u';
            nu = obj.c(3);
            AA = 2*(2*nu-1)/(nu-1);
            penalty = 1e10;
            N = size(x,1);
            A = zeros(2*N,2*N); % 2N GDLs
            m = weightedVolume(par_omega)*ones(length(x),1);
            for ii = 1:N
                c1_hat = (obj.c(1) + obj.c(2)/9*m(ii)*(-nu+2)/(2*nu-1));
                dofi = [idb(2*ii-1) idb(2*ii)];
                family = familySet(ii,familySet(ii,:)~=0);
                Nj = length(family);
                iII = 1:Nj;
                jj = family'; % j sum
                dofj = [idb(2*jj-1) idb(2*jj)];
                eta = u(dofj) - u(dofi);
                xi = x(jj,:) - x(ii,:);
                M = eta+xi;
                normaj = vecnorm(xi,2,2); 
                norma_eta = vecnorm(M,2,2);
                omegaj = influenceFunction(normaj,par_omega);
                muj = mu{ii};
                PSI_ij = M./norma_eta;
                % Parameters
                Vij = partialAreas(ii,iII)';
                c1 = obj.c(2)*omegaj;
                c2 = c1_hat*AA/m(ii).*omegaj.*normaj;
                g = omegaj.*normaj;
                % U and V
                if dofi(1) <= ndof || dofi(2) <= ndof
                    % First dof of node ii is free
  
                    ti1 = sum(-2*c1.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij) + sum(AA/m(ii).*g.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj).*sum(-c2.*Vij.*PSI_ij); % Aii
                    ti2 = sum(-2*c1.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij) + sum(AA/m(ii).*g.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj).*sum(-c2.*Vij.*PSI_ij); % Aip
                    tj1 = 2*c1.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij + AA/m(ii).*(g.*PSI_ij(:,1).*Vij.*surfaceCorrection(ii,iII)'.*muj)*sum(c2.*Vij.*PSI_ij);% Aij
                    tj2 = 2*c1.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj.*PSI_ij + AA/m(ii).*(g.*PSI_ij(:,2).*Vij.*surfaceCorrection(ii,iII)'.*muj)*sum(c2.*Vij.*PSI_ij);% Aijp
                    
                    for Ij = 1:length(jj)
                        j = jj(Ij);
                        kk  = familySet(j,familySet(j,:)~=0)';
                        iIII = 1:length(kk);
                        dofk = [idb(2*kk-1) idb(2*kk)];
                        eta_k = u(dofk) - u(dofj(Ij,:));
                        xi_k = x(kk,:) - x(j,:);
                        M_k = eta_k+xi_k;
                        normak = vecnorm(xi_k,2,2); 
                        
                        norma_etak = vecnorm(M_k,2,2);
                        omegak = influenceFunction(normak,par_omega);
                        muk = mu{j};
                        PSI_jk = M_k./norma_etak;
                        Vjk = partialAreas(j,iIII)';
                        % Parameters
                        gk = omegak.*normak;
                        
                        tj1(Ij,:) = tj1(Ij,:) - c2(Ij).*sum(AA/m(j).*gk.*PSI_jk(:,1).*Vjk.*surfaceCorrection(j,iIII)'.*muk).*Vij(Ij).*PSI_ij(Ij,:); 
                        tj2(Ij,:) = tj2(Ij,:) - c2(Ij).*sum(AA/m(j).*gk.*PSI_jk(:,2).*Vjk.*surfaceCorrection(j,iIII)'.*muk).*Vij(Ij).*PSI_ij(Ij,:);
                        tk1 = c2(Ij)*AA/m(j).*(gk.*PSI_jk(:,1).*Vjk.*surfaceCorrection(j,iIII)'.*muk)*Vij(Ij).*PSI_ij(Ij,:);
                        tk2 = c2(Ij)*AA/m(j).*(gk.*PSI_jk(:,2).*Vjk.*surfaceCorrection(j,iIII)'.*muk)*Vij(Ij).*PSI_ij(Ij,:);
                        
                        if dofi(1) <= ndof
                            A(dofi(1),dofk(:,1)) = A(dofi(1),dofk(:,1)) + tk1(:,1)' * V(ii);
                            A(dofi(1),dofk(:,2)) = A(dofi(1),dofk(:,2)) + tk2(:,1)' * V(ii);
                        end
                        if dofi(2) <= ndof
                            A(dofi(2),dofk(:,1)) = A(dofi(2),dofk(:,1)) + tk1(:,2)' * V(ii);
                            A(dofi(2),dofk(:,2)) = A(dofi(2),dofk(:,2)) + tk2(:,2)' * V(ii);
                        end 
                    end
                    % U
                    if dofi(1) <= ndof
                        A(dofi(1),dofi(1)) = A(dofi(1),dofi(1)) + ti1(1) *V(ii) ;
                        A(dofi(1),dofj(:,1)) = A(dofi(1),dofj(:,1)) + tj1(:,1)' * V(ii);
                        A(dofi(1),dofi(2)) = A(dofi(1),dofi(2)) + ti2(1) * V(ii);
                        A(dofi(1),dofj(:,2)) = A(dofi(1),dofj(:,2)) + tj2(:,1)'* V(ii);
                    else
                        % Constraint nodes
                        A(dofi(1),dofi(1)) = -penalty;
                    end
                   % V
                   if dofi(2) <= ndof
                        A(dofi(2),dofi(1)) = A(dofi(2),dofi(1)) + ti1(2) *V(ii) ;
                        A(dofi(2),dofj(:,1)) = A(dofi(2),dofj(:,1)) + tj1(:,2)' * V(ii);
                        A(dofi(2),dofi(2)) = A(dofi(2),dofi(2)) + ti2(2) * V(ii);
                        A(dofi(2),dofj(:,2)) = A(dofi(2),dofj(:,2)) + tj2(:,2)' * V(ii);
                   else
                        % Constraint nodes
                        A(dofi(2),dofi(2)) = -penalty;
                   end
                else
                    A(dofi(1),dofi(1)) = -penalty;
                    A(dofi(2),dofi(2)) = -penalty;
                end
            end
        end
    end
end
function ff = fscalar(obj,x,damage,noFail,ii)

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

    % NEW FORMULATION
    ff = (x>S1(1)).*(x<S0(1)).*(S0(1)*(x-S1(1))./(S0(1) - S1(1))) + (x>=S0(1)).*(x<=S0(2)).*x ...
        + (x>S0(2)).*(x<S1(2)).*(S0(2)*(S1(2)-x)./(S1(2)-S0(2)));
    % Adding noFail condition
    ff(noFail) = x(noFail);
else
    ff = x;
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

