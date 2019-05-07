function [f,history] = interactionForce_LLPS(x,u,ii,dof_vec,familyMat,partialArea,neighIndex,history,dt,noFail)
    global c1 c2 horizon omega Sc 
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
    m = weightedVolume(horizon,omega);
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
            theta(bond == ii + 2*bond == jj)= 2/m*influenceFunction(norm(zeta),horizon,omega)*S_z*norm(zeta)^2*partialArea(bond,neighInd2) +  theta(bond == ii + 2*bond == jj); % The boolean operator returns 1 if equal to ii and 2 if equal to jj
            neighInd2 = neighInd2 + 1;
        end
    end
    T_ii = c1*influenceFunction(norma,horizon,omega)*theta(1)*xi + c2*influenceFunction(norma,horizon,omega)*ext*ee;
    T_jj = c1*influenceFunction(norma,horizon,omega)*theta(2)*(-xi) + c2*influenceFunction(norma,horizon,omega)*(-ext)*ee;
    f = T_ii - T_jj;
end