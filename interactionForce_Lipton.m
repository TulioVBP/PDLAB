function [f,history_up,mu] = interactionForce_Lipton(x,u,ii,dof_vec,familyMat,partialAreas,neighIndex,separatorDamage,dt,history,noFail)
%% INPUT
% x_i: position of node i
% x_j: position of node j
% u_i: displacement of node i
% u_j: displacement of node j
% dt: time step
% S_max_ant: maximum stretch for each given bond
% notch: coordinates of the initial notch
% noFail: true if the damage is off for this specific bond
%% OUTPUT
% f: vector internal force between j and i nodes
% S_max: maximum stretch for each bond
% T: vector state force
%% CODE
    global horizon omega Sc
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
    js = @(x) (x >= Sc).*(x/Sc-1).^5./(1+(x/Sc).^5); % If x < Sc, then js = 0;
    history_up(1) = history(1) + js(S)*dt;
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
        theta_i = theta_i + 1/V_delta*influenceFunction(norm(zeta),horizon,omega)*norm(zeta)^2*S_z*partialAreas(ii,neighIndex2);
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
        theta_j = theta_j + 1/V_delta*influenceFunction(norm(zeta),horizon,omega)*norm(zeta)^2*S_z*partialAreas(jj,neighIndex2);
        neighIndex2 = neighIndex2  + 1;
    end
    thetac_p = 0.01;
    thetac_m = 0.01;
    jth = @(x) (x >= thetac_p).*(x/thetac_p-1).^4./(1+(x/thetac_p).^5) + (x <= -thetac_m).*(-x/thetac_m-1).^4./(1+(-x/thetac_m).^5); % If x < Sc, then js = 0;
    history_up(2) = history(2) + jth(theta_i)*dt; % Update integral of dilatation of x_i for this specific interaction
    history_up(3) = history(3) + jth(theta_j)*dt; % Update integral of dilatation of x_j for this specific interaction
    H = damageFactor(history_up,x_i,x_j,noFail);
    Ht = H(1); Hd_x = H(2); Hd_y = H(3);
    % Tensile term Lt
    ft = 2/V_delta*influenceFunction(norma,horizon,omega)/horizon*Ht*fscalar(sqrt(norma)*S,norma)*ee;
    % Dilatation term Ld
    fd = 1/V_delta*influenceFunction(norma,horizon,omega)/horizon^2*norma*(Hd_y*gscalar(theta_j) + Hd_x*gscalar(theta_i))*ee;
    % Final force
    f = fd + ft;
    mu = Ht; % Check for damage in this model
end

function ff = fscalar(x,norma)
global damageOn c1
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x <= r1
        ff = c1*x*sqrt(norma);
    elseif x > r2
        ff = norma;
    end
else
    ff = c1*x*sqrt(norma);
end
end

function gg = gscalar(x)
global damageOn c2
r1 = 3.0;
r2 = 3.0;
if damageOn
    if x <= r1
        gg = c2*x;
    elseif x > r2
        gg = 1;
    end
else
    gg = c2*x;
end

end
