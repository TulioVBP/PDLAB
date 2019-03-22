function phi = damageIndex(x,u,family,partialAreas,ii)
    x_i = x(ii,:);
    u_i = u(ii,:);
    partialDamage = 0; % Total damage
    areaTot = 0; % Total area of the neighbourhood
    for familyIndex = 1:length(family)
       if family(familyIndex)~= 0
           jj = family(familyIndex);
           x_j = x(jj,:);
           u_j = u(jj,:);
           xi = x_j - x_i; % \xi
           eta = u_j - u_i; % \eta
           norma = norm(xi); 
           S = (norm(eta+xi) - norm(xi))/norma; % Calculate stretch
           mu_j = damageFactor(S); % Calculate damage factor
           areaTot = areaTot + partialAreas(familyIndex); % Integration on the area
           partialDamage = partialDamage + mu_j*partialAreas(familyIndex); % Integration on the partial damage
       else
           break;
       end
    end
    phi = 1 - partialDamage/areaTot;
end