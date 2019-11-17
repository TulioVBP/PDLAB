function [f,history,mu] = interactionForce_LLPS(x,u,ii,dof_vec,familyMat,partialArea,neighIndex,par_omega,c,model,separatorDamage,damage,dt,history,noFail) 
    jj = familyMat(ii,neighIndex);
    x_i = x(ii,:); x_j = x(jj,:);
    dofi = dof_vec(ii,:); dofj = dof_vec(jj,:);
    xi = x(jj,:) - x(ii,:); % \xi
    u_i = u(dofi)'; u_j = u(dofj)';
    eta = u_j - u_i; % \eta
    norma = norm(xi); 
    ext = dot(eta,xi)/norma; % Calculate stretch - linear
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
    % ---- Evaluatin the force ----
    % Dilatation term
    theta = [0 0];
    m = weightedVolume(par_omega);
    for bond = [ii jj]
        dofb = dof_vec(bond,:);
        u_b = u(dofb)';
        neighInd2 = 1;
        for kk = familyMat(bond,familyMat(bond,:)~=0)
            zeta = x(kk,:) - x(bond,:);
            dofk = dof_vec(kk,:);
            u_k = u(dofk)';
            eta_2 = u_k - u_b;
            S_z = dot(zeta,eta_2)/norm(zeta)^2;
            theta(bond == ii + 2*bond == jj)= 2/m*influenceFunction(norm(zeta),par_omega)*S_z*norm(zeta)^2*partialArea(bond,neighInd2) +  theta(bond == ii + 2*bond == jj); % The boolean operator returns 1 if equal to ii and 2 if equal to jj
            neighInd2 = neighInd2 + 1;
        end
    end
    T_ii = c(1)*influenceFunction(norma,par_omega)*theta(1)*xi + c(2)*influenceFunction(norma,par_omega)*ext*ee;
    T_jj = c(1)*influenceFunction(norma,par_omega)*theta(2)*(-xi) + c(2)*influenceFunction(norma,par_omega)*(-ext)*ee;
    f = T_ii - T_jj;
    mu = 1; % No damage in this model
end