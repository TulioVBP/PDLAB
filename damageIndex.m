function phi = damageIndex(x,u,family,partialAreas,ii,notch,idb,noFailZone)
%% INPUT:
% - x: nodes position
% - u: nodes displacement
% - family: family of node ii
% - partialAreas: partial areas of nodes jj inside ii neighborhood
% - ii: i-th node used
% - notch: initial notch present in the material
% - idb: indexing each dof to its position
% - noFailZone: array with zeros and ones
%% OUTPUT:
% phi: damage index for the i-th node
%% CODE
    x_i = x(ii,:);
    oneDcond = length(u) == 2*length(x);
    if oneDcond
        dofi = [idb(2*ii-1) idb(2*ii)];
        u_i = u(dofi)';
    else
        u_i = u(ii,:);
    end
    partialDamage = 0; % Total damage
    areaTot = 0; % Total area of the neighbourhood
    for familyIndex = 1:length(family)
       if family(familyIndex)~= 0
           jj = family(familyIndex);
           x_j = x(jj,:);
           if oneDcond
            dofj = [idb(2*jj-1) idb(2*jj)];
            u_j = u(dofj)';
           else
            u_j = u(jj,:);
           end
           xi = x_j - x_i; % \xi
           eta = u_j - u_i; % \eta
           norma = norm(xi); 
           S = (norm(eta+xi) - norm(xi))/norma; % Calculate stretch
           noFail = noFailZone(ii) || noFailZone(jj);
           mu_j = damageFactor(S,notch,x_i,x_j,noFail); % Calculate damage factor
           areaTot = areaTot + partialAreas(familyIndex); % Integration on the area
           partialDamage = partialDamage + mu_j*partialAreas(familyIndex); % Integration on the partial damage
       else
           break;
       end
    end
    phi = 1 - partialDamage/areaTot;
end