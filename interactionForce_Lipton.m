function [f,history_up,mu] = interactionForce_Lipton(x,u,ii,dof_vec,familyMat,partialAreas,neighIndex,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
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
    jj = familyMat(ii,neighIndex);
    x_i = x(ii,:); x_j = x(jj,:);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x(jj,:) - x(ii,:); % \xi
    u_i = u(dofi)'; u_j = u(dofj)';
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    S = dot(eta,xi)/norma^2; % Calculate stretch - linear
    ee = (xi)/norma; % Versor
    % Defining non-existing parameters
    if ~exist('noFail','var')
        noFail = 0;
    end
    if ~exist('history','var')
        history = [0 0 0];
    end
    if ~exist('dt','var')
       dt = 0; 
    end
    %% Evaluating the force interaction
    V_delta = pi*horizon^2; % Not sure if this is the right expression
    % ---- Evaluating the damage factor Ht
    % . Evaluating js
    %js = @(x) (x >= Sc).*(x/Sc-1).^5./(1+(x/Sc).^5); % If x < Sc, then js = 0;
    history_up(1) = history(1) + js(S,damage.Sc)*dt;
    % -- Evaluating the dilatation
    % Dilatation term
    theta_i = 0;
    neighIndex2 = 1;
    for kk = familyMat(ii,familyMat(ii,:)~=0)
        zeta = x(kk,:) - x(ii,:);
        dofk = dof_vec(kk,:);
        u_k = u(dofk)';
        eta_2 = u_k - u_i;
        S_z = dot(zeta,eta_2)/norm(zeta)^2;
        theta_i = theta_i + 1/V_delta*influenceFunction(norm(zeta),par_omega)*norm(zeta)^2*S_z*partialAreas(ii,neighIndex2);
        neighIndex2 = neighIndex2  + 1;
    end
    theta_j = 0;
    neighIndex2 = 1;
    for kk = familyMat(jj,familyMat(jj,:)~=0)
        zeta = x(kk,:) - x(jj,:);
        dofk = dof_vec(kk,:);
        u_k = u(dofk)';
        eta_2 = u_k - u_j;
        S_z = dot(zeta,eta_2)/norm(zeta)^2;
        theta_j = theta_j + 1/V_delta*influenceFunction(norm(zeta),par_omega)*norm(zeta)^2*S_z*partialAreas(jj,neighIndex2);
        neighIndex2 = neighIndex2  + 1;
    end
    thetac_p = 0.01;
    thetac_m = 0.01;
    %jth = @(x) (x >= thetac_p).*(x/thetac_p-1).^4./(1+(x/thetac_p).^5) + (x <= -thetac_m).*(-x/thetac_m-1).^4./(1+(-x/thetac_m).^5); % If x < Sc, then js = 0;
    history_up(2) = history(2) + jth(theta_i,thetac_p,thetac_m)*dt; % Update integral of dilatation of x_i for this specific interaction
    history_up(3) = history(3) + jth(theta_j,thetac_p,thetac_m)*dt; % Update integral of dilatation of x_j for this specific interaction
    H = damageFactor(history_up,x_i,x_j,damage,noFail,model);
    Ht = H(1); Hd_x = H(2); Hd_y = H(3);
    % Tensile term Lt
    ft = 2/V_delta*influenceFunction(norma,par_omega)/horizon*Ht*fscalar(sqrt(norma)*S,norma,c,damage.damageOn)*ee;
    % Dilatation term Ld
    fd = 1/V_delta*influenceFunction(norma,par_omega)/horizon^2*norma*(Hd_y*gscalar(theta_j,c,damage.damageOn) + Hd_x*gscalar(theta_i,c,damage.damageOn))*ee;
    % Final force
    f = fd + ft;
    mu = Ht; % Check for damage in this model
end

function ff = fscalar(x,norma,c,damageOn)
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x*sqrt(norma) <= r1
        ff = c(1)*x*sqrt(norma);
    elseif x*sqrt(norma) > r2
        ff = sqrt(norma);
    end
else
    ff = c(1)*x*sqrt(norma);
end
end

function gg = gscalar(x,c,damageOn)
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x <= r1
        gg = c(2)*x;
    elseif x > r2
        gg = 1;
    end
else
    gg = c(2)*x;
end

end

function jj = js(x,Sc)
    jj = (x >= Sc).*(x/Sc-1).^5./(1+(x/Sc).^5); % If x < Sc, then js = 0;
end

function jj = jth(x,thetac_p,thetac_m)
    jj = (x >= thetac_p).*(x/thetac_p-1).^4./(1+(x/thetac_p).^5) + (x <= -thetac_m).*(-x/thetac_m-1).^4./(1+(-x/thetac_m).^5); % If x < Sc, then js = 0;
end