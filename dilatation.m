function [theta, history_thetaUp] = dilatation(x,u,family,partialAreas,surfaceCorrection,transvList,idb,par_omega,c,model,damage,history_theta,dt)
    horizon = par_omega(1);
    thetac_p = 0.01; % For Lipton's damage
    thetac_m = 0.01; % For Lipton's damage
    m = weightedVolume(par_omega);
    if isempty(transvList) % Not a specific range of nodes was chosen
       transvList = 1:length(x);
    end
    if exist('history_theta','var')
        history_thetaUp = history_theta;
    else
        history_thetaUp = [];
    end
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
        norma = vecnorm(xi')';
        eta = u(dofj)-u(dofi);
            switch model.number            
                case 6%"Linearized LPS"
                    theta_vec = 3/m*influenceFunction(norma,par_omega).*dot(eta',xi')'.*partialAreas(transv_ind,neigh_ind)'.*surfaceCorrection(transv_ind,neigh_ind)';
                    theta(transv_ind) = sum(theta_vec);
                case 3 %"Lipton Free Damage"
                    V_delta = pi*horizon^2;
                    S_linear = dot(xi',eta')'./norma.^2;
                    theta_vec = 1/V_delta*influenceFunction(norma,par_omega).*norma.^2.*S_linear.*partialAreas(transv_ind,neigh_ind)'.*surfaceCorrection(transv_ind,neigh_ind)';
                    if nargin > 10
                        wholeBonds = ~damage.brokenBonds(ii,neigh_ind)';
                        theta(transv_ind) = sum(theta_vec.*wholeBonds);
                    else
                        theta(transv_ind) = sum(theta_vec);
                    end
                case 4 %"LPS 2D"
                    nu = c(3);
                    elong = vecnorm(xi'+eta')' - vecnorm(xi')';
                    theta_vec = 2*(2*nu-1)/(nu-1)/m*influenceFunction(norma,par_omega).*norma.*elong.*partialAreas(transv_ind,neigh_ind)'.*surfaceCorrection(transv_ind,neigh_ind)';
                    theta(transv_ind) = sum(theta_vec);    
                otherwise
                    break;
            end
        if nargin == 12
            history_thetaUp(transv_ind) = history_theta(transv_ind) + jth(theta(transv_ind),thetac_p,thetac_m)*dt; % Update integral of dilatation of x_i for this specific interaction 
        end
        transv_ind = transv_ind + 1; 
    end  
end

function jj = jth(x,thetac_p,thetac_m)
    jj = (x >= thetac_p).*(x/thetac_p-1).^4./(1+(x/thetac_p).^5) + (x <= -thetac_m).*(-x/thetac_m-1).^4./(1+(-x/thetac_m).^5); % If x < Sc, then js = 0;
end
