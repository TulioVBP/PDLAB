function [f,history,mu] = interactionForce_StateBased(x,u,ii,dof_vec,familyMat,partialAreas,neighIndex,par_omega,c,model,separatorDamage,damage,dt,history,noFail)
    nu = c(3);
    jj = familyMat(ii,neighIndex);
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
            e_z = norm(zeta+eta_2) - norm(zeta); % Elongation
            theta((bond == ii) + 2*(bond == jj))= 2*(2*nu-1)/(nu-1)/m*influenceFunction(norm(zeta),par_omega)*norm(zeta)*e_z*partialAreas(bond,neighInd2) +  theta((bond == ii) + 2*(bond == jj)); % The boolean operator returns 1 if equal to ii and 2 if equal to jj
            neighInd2 = neighInd2 + 1;
        end
    end
    % Additional integral term
    add_term = [0 0];
    for bond = [ii jj]
        dofb = dof_vec(bond,:);
        u_b = u(dofb)';
        neighInd2 = 1;
        for kk = familyMat(bond,familyMat(bond,:)~=0)
            zeta = x(kk,:) - x(bond,:);
            dofk = dof_vec(kk,:);
            u_k = u(dofk)';
            eta_2 = u_k - u_b;
            e_z = norm(zeta+eta_2) - norm(zeta); % Elongation
            ed = e_z - theta((bond == ii) + 2*(bond == jj))*norm(zeta)/3;
            add_term((bond == ii) + 2*(bond == jj)) = influenceFunction(norm(zeta),par_omega)*norm(zeta)*ed*partialAreas(bond,neighInd2) + add_term((bond == ii) + 2*(bond == jj)); 
            neighInd2 = neighInd2 + 1;
        end
    end
    % State forces
    T_ii = 2*(2*nu-1)/(nu-1)/m*(c(1)*theta(1) - c(2)/3*add_term(1))*influenceFunction(norma,par_omega)*norma ...
        + c(2)*influenceFunction(norma,par_omega)*(elong - theta(1)*norma/3);%c1*influenceFunction(norma,horizon,omega)*theta(1)*xi + c2*influenceFunction(norma,horizon,omega)*ext*ee;
    T_jj = 2*(2*nu-1)/(nu-1)/m*(c(1)*theta(2) - c(2)/3*add_term(2))*influenceFunction(norma,par_omega)*norma ...
        + c(2)*influenceFunction(norma,par_omega)*(elong - theta(2)*norma/3);
    f = (T_ii + T_jj)*ee;
    mu = 1; % No damage in this model
end