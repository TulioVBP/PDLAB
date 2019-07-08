function [f,historyS,mu] = interactionForce_StateBased(x,u,theta,ii,jj,dof_vec,par_omega,c,model,separatorDamage,damage,dt,historyS,historyTheta,noFail)
    nu = c(3);
    x_i = x(ii,:); x_j = x(jj,:);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x(jj,:) - x(ii,:); % \xi
    u_i = u(dofi)'; u_j = u(dofj)';
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    elong = norm(xi+eta) - norma; % Calculate elongation - linear
    ee = (xi+eta)/norm(xi+eta); % Versor
    % Defining non-existing parameters
    if ~exist('noFail','var')
        noFail = 0;
    end
    if ~exist('history','var')
        historyS = 0;
        historyTheta = zeros(length(x),1);
    end
    if ~exist('dt','var')
       dt = 0; 
    end
    % ---- Evaluatin the force ----
    % Dilatation term
    theta = [theta(ii) theta(jj)];
    m = weightedVolume(par_omega);
    T_ij = 2*(2*nu-1)/(nu-1)*((c(1) + c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega)*norma/m)*theta(1) ...
        + c(2)*influenceFunction(norma,par_omega)*(elong);
    T_ji = 2*(2*nu-1)/(nu-1)*((c(1) + c(2)/9*m*(-nu+2)/(2*nu-1))*influenceFunction(norma,par_omega)*norma/m)*theta(2) ...
        + c(2)*influenceFunction(norma,par_omega)*(elong);
    f = (T_ij + T_ji)*ee;
    mu = 1; % No damage in this model
    historyS = 0;
end