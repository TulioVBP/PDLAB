function [f,history_upS,mu] = interactionForce_LSJT(x,u,theta,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,historyS,historyTheta,noFail)
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
    %jj = familyMat(ii,neighIndex);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x(jj,:) - x(ii,:); % \xi
    u = u';
    u_i = u(dofi); u_j = u(dofj);
    eta = u_j - u_i; % \eta
    norma = vecnorm(xi')'; 
    S = dot(eta',xi')'./norma.^2; % Calculate stretch - linear
    ee = (xi)./norma; % Versor
    if nargin > 10  && damage.damageOn% Damage considered
        % ---- Evaluating the damage factor Ht
        % . Evaluating js
        history_upS = historyS' + js(S,damage.Sc)*dt; % NX1
        XX = [history_upS, historyTheta(ii)*ones(length(historyTheta(jj)),1), historyTheta(jj)]; 
        H = damageFactor(XX,ii,1:length(jj),damage,noFail,model);
    else
        history_upS = 0;
        H = ones(length(jj),3);
    end
    %% Evaluating the force interaction
    V_delta = pi*horizon^2; % Not sure if this is the right expression
    
    % -- Evaluating the dilatation
    % Dilatation term
    theta_i = theta(ii); theta_j = theta(jj);
    
    Ht = H(:,1); Hd_x = H(1,2); Hd_y = H(1,3);
    % Tensile term Lt
    ft = 2/V_delta*influenceFunction(norma,par_omega)./horizon.*Ht.*fscalar(sqrt(norma).*S,norma,c,damage.damageOn,damage.Sc).*ee;
    % Dilatation term Ld
    fd = 1/V_delta*influenceFunction(norma,par_omega)./horizon^2.*norma.*(Hd_y.*gscalar(theta_j,c,damage.damageOn,damage.thetaC) + Hd_x.*gscalar(theta_i,c,damage.damageOn,damage.thetaC)).*ee;
    % Final force
    f = fd + ft;
    mu = Ht; % Check for damage in this model
    if any(damage.brokenBonds(ii,1:length(jj)))
        f(damage.brokenBonds(ii,1:length(jj)),:) = zeros(sum(damage.brokenBonds(ii,1:length(jj))),2); % Guaranteeing initial damage
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
    gg = (x<=r1).*c(2).*x + (x>2)*1;
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