classdef modelLSJT
     properties
        linearity = true;
        stiffnessAnal = false;
        dilatation = true;
        damage;
        c;
        history;
        theta;
        dt;
     end
    methods
        function obj = modelLSJT(E,nu,par_omega,damage,G0,dt)
            %% Pre Initialization %%
             % Any code not using output argument (obj)
             if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                E = 1;
                nu = 1/4;
                par_omega = [1 1 1];
                dt = 1;
                damage.DD = false;
                G0 = 0;
             end
            horizon = par_omega(1);
            mm = weightedVolume(par_omega);
            obj.c(1) = 8*pi*horizon^3/mm*E/(1+nu)/2;
            obj.c(2) = -2*(-1 + 3*nu)/(-1 + nu^2)*(pi*horizon^3/mm)^2*E;
            obj.dt = dt;
            if par_omega(2) == 3 && par_omega(3) == 1
                obj.damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E)/horizon);
            elseif par_omega(2) == 1 && par_omega(3) == 0
                l = 3;
                obj.damage.Sc = sqrt((1+1/3)*pi*G0*l/(8*E*horizon*0.66467));
            elseif par_omega(2) == 1 && par_omega(3) == 1
                l = 3;
                obj.damage.Sc = sqrt(G0*(1/3+1)*pi^(3/2)/8/E*(l/horizon));
            else
                warning('Critical bond not defined.')
            end
            obj.damage.thetaC_p = 3; % Very high value to not be activated
            obj.damage.thetaC_m = -3; % Very low value to not be activated 
            % Damage dependent Sc
            if damage.DD
                obj.damage.alfa = 0.2; obj.damage.beta = 0.2; obj.damage.gamma = 1.4;
            else
                obj.damage.alfa = 0; obj.damage.beta = 0; obj.damage.gamma = 1; % No dependency
            end
        end
        
        function [theta,historyT] = dilatationEval(obj,x,u,family,partialAreas,surfaceCorrection,transvList,idb,par_omega,damage,historyS,historyT)
            if isempty(transvList) % Not a specific range of nodes was chosen
               transvList = 1:length(x);
            end
            horizon = par_omega(1);
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
                V_delta = pi*horizon^2;
                S_linear = dot(xi',eta')'./norma.^2;
                theta_vec = 1/V_delta*influenceFunction(norma,par_omega).*norma.^2.*S_linear.*partialAreas(transv_ind,neigh_ind)'.*surfaceCorrection(transv_ind,neigh_ind)';
                if damage.damageOn
                    Sc = obj.damage_dependent(damage.phi(ii));
                    history_upS = obj.updateHistory(S_linear,historyS(ii,neigh_ind)',Sc);
                    XX = {history_upS, historyT(ii), historyT(jj)}; 
                    noFail = damage.noFail(ii) | damage.noFail(jj); % True if node ii or jj is in the no fail zone
                    [H,~,~] = obj.damageFactor(XX,ii,1:length(jj),damage,noFail);
                    theta(transv_ind) = sum(theta_vec.*H); % Tulio's model
                end
                transv_ind = transv_ind + 1;
            end
            historyT = obj.updateHistoryT(theta(transvList),historyT);
        end
        
        function [f,historyS,mu] = T(obj,x,u,theta,ii,jj,dof_vec,par_omega,separatorDamage,damage,historyS,historyTheta,noFail)
            %% INPUT
        % x: nodes position
        % u: dof's displacement
        % ii: i-th node index
        % dof_vec: index of k-th(k is the row index) dof's
        % familyMat: family interaction matrix
        % partialAreas: partial areas 
        % neighIndex: index of the j-th node
        % par_omega: parameters for the omega
        % c: material's constants
        % separatorDamage: doesn't have a role but it helps differentiate between
        %                  what is damage related and what is not
        % dt: step time 
        % history: depending on the model, is an arbitrary history dependent
        %          variable
        % noFail: true if the either of the nodes jj and ii is a no-fail zone

        %% OUTPUT
        % f: vector internal force between j and i nodes
        % history_up: updated history dependent variable
        % T: vector state force
        %% CODE
            horizon = par_omega(1);
            dofi = dof_vec(ii,:); 
            dofj = dof_vec(jj,:);
            xi = x(jj,:) - x(ii,:); % \xi
            u = u';
            u_i = u(dofi); 
            u_j = u(dofj);
            eta = u_j - u_i; % \eta
            norma = vecnorm(xi')'; 
            S = dot(eta',xi')'./norma.^2; % Calculate stretch - linear
            ee = (xi)./norma; % Versor
            if nargin > 10  && damage.damageOn% Damage considered
                % . Evaluating js
                Sc = obj.damage_dependent(damage.phi(ii));
                historyS = obj.updateHistory(S,historyS',Sc); % NX1
                XX = {historyS, historyTheta(ii), historyTheta(jj)}; 
                [Ht,Hd_x,Hd_y] = obj.damageFactor(XX,ii,1:length(jj),damage,noFail);
            else
                historyS = 0;
                H = ones(length(jj),3);
            end
            historyS = historyS';
            %% Evaluating the force interaction
            V_delta = pi*horizon^2; % Not sure if this is the right expression

            % -- Evaluating the dilatation
            % Dilatation term
            theta_i = theta(ii); theta_j = theta(jj);


            % Tensile term Lt
            ft = 2/V_delta*influenceFunction(norma,par_omega)./horizon.*Ht.*fscalar(sqrt(norma).*S,norma,obj.c,damage.damageOn,Sc).*ee;
            % Dilatation term Ld
            fd = 1/V_delta*influenceFunction(norma,par_omega)./horizon^2.*norma.*(Hd_y.*gscalar(theta_j,obj.c,damage.damageOn,obj.damage.thetaC_p) + Hd_x.*gscalar(theta_i,obj.c,damage.damageOn,obj.damage.thetaC_p)).*ee;
            % Final force
            f = fd + ft;
            mu = Ht; % Check for damage in this model
            if any(damage.brokenBonds(ii,1:length(jj)))
                f(damage.brokenBonds(ii,1:length(jj)),:) = zeros(sum(damage.brokenBonds(ii,1:length(jj))),2); % Guaranteeing initial damage
            end
        end
        
        function [historyS] = updateHistory(obj,S,historyS,Sc)
            if nargin < 2
                error('No history was passed as argument.')
            end
            % Updating history values
            historyS = historyS + js(S,Sc)*obj.dt;
        end
        
        function [history_theta] = updateHistoryT(obj,theta,history_theta)
             history_theta = history_theta + jth(theta,obj.damage.thetaC_p,obj.damage.thetaC_m)*obj.dt;
        end
        
        function [HT,HDi,HDj] = damageFactor(obj,x,ii,neighIndex,damage,noFail)
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
            
            % HT
            xc = (0.05)^2/(1+1.05^2) * 0.02e-6; % js(Sc)*dt = 2.3781e-11
            HT = (x{1}<xc).*(exp(1-1./(1-(x{1}/xc).^2.01)));
            HT(isnan(HT)) = zeros(sum(sum(isnan(HT))),1);
            % {HD ALWAYS ONE}
            HDi = 1;
            HDj = ones(length(x{3}),1);
            % Deal with initial broken bonds
            brokenBonds = damage.brokenBonds(ii,neighIndex);
            if any(brokenBonds)
                HT(brokenBonds,:) = zeros(sum(brokenBonds),size(HT,2));
            end
            if ~isempty(noFail)
                HT(noFail,:) = ones(size(HT(noFail,:)));
            end
        end
        
        function W = strainEnergyDensity(obj,x,u,theta,family,partialAreas,surfaceCorrection,ii,idb,par_omega,damage,historyS,historyT) 
            familySet = family(family~=0);
            dofi = [idb(2*ii-1) idb(2*ii)];
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
            if nargin > 11 && damage.damageOn
                noFail = damage.noFail(ii) | damage.noFail(jj);
                XX = {historyS(neigh_ind)', historyT(ii)*ones(length(jj),1),historyT(jj)};
                [Ht,Hd_x,~] = obj.damageFactor(XX,ii,neigh_ind,damage,noFail);
                Sc = obj.damage_dependent(damage.phi(ii));
            else
                Ht = 1;
                Hd_x = 1;
                Sc = obj.damage.Sc;
            end
            Slin = dot(xi',eta')'./norma.^2;
            V_delta = pi*horizon^2;
            w = 1/V_delta*(influenceFunction(norma,par_omega).*norma/horizon.*Ht.*f_potential(Slin,sqrt(norma),obj.c,damage,Sc));
            W = sum(w.*partialAreas(neigh_ind)'.*surfaceCorrection(neigh_ind)') + 1/horizon^2* Hd_x *g_potential(theta_i,obj.c,damage,obj.damage.thetaC_m);
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
                obj.history.S = values{1};
            end
        end
        
        function Sc = damage_dependent(obj,phi)
            %Damage dependent crack
            alfa = obj.damage.alfa; beta = obj.damage.beta; gamma = obj.damage.gamma;
            if phi > alfa
                Sc = obj.damage.Sc*min(gamma,1+beta*(phi-alfa)/(1-phi));
            else
                Sc = obj.damage.Sc;
            end
        end
    end
end

function ff = fscalar(x,norma,c,damageOn,Sc)
r1 = Sc*sqrt(norma);
%r1 = 3; % Uncomment for a better result
r2 = r1;
if damageOn
    ff = (x.*sqrt(norma) <= r1).*c(1).*x.*sqrt(norma)...
        + (x.*sqrt(norma) > r2).*sqrt(norma);
else
    ff = c(1)*x.*sqrt(norma);
end
end

function gg = gscalar(x,c,damageOn,thetaC)
r1 = thetaC;
r2 = r1;
if damageOn
    %gg = (x<=r1).*c(2).*x + (x>2)*1;
    gg = c(2)*x;
else
    gg = c(2)*x;
end

end

function jj = js(x,Sc)
    jj = (x >= Sc).*(x/Sc-1).^2./(1+(x/Sc).^2); % If x < Sc, then js = 0;
end

function jj = jth(x,thetac_p,thetac_m)
    jj = (x >= thetac_p).*(x/thetac_p-1).^4./(1+(x/thetac_p).^5) + (x <= -thetac_m).*(-x/thetac_m-1).^4./(1+(-x/thetac_m).^5); % If x < Sc, then js = 0;
end

function ff = f_potential(S,norma,c,damage,Sc)
r1 = Sc.*norma;
%r1 = 3; % Uncomment for a better result
r2 = r1;
x = S.*norma;
if damage.damageOn
    ff = (x <= r1).*c(1).*x.^2/2 + (x > r2).*x;
else
    ff = c(1)*x.^2/2;
end
end

function gg = g_potential(x,c,damage,thetaC)
r1 = thetaC;
r2 = r1;
if damage.damageOn
    %gg = (x <= r1).*c(2).*x.^2/2 + (x > r2).*x;
    gg = c(2).*x.^2/2;
else
    gg = c(2)*x^2/2;
end

end

