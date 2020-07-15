function [model,c,T,damage] = modelParameters(model,par_omega,damage,E,nu,G0)
% INPUTS
% - model: the constitutive model you are using. Implemented models are:
% "DTT", "PMB", "LBB", "LSJ", "LSJ-T","LPS-T"
% - par_omega = [horizon omega alfa]
% - damage: struct variable with damage parameters
% - E: Young modulus
% - nu: Poisson's ratio
% - G0: Energy release rate
switch model.name
    case "DTT"
        %% PMB DTT
        alfa = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(par_omega);
        horizon = par_omega(1);
        c(1) = 6*E/mm;
        T = @models.forces.interactionForce_PMBDTT;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
        % --- SC
        if par_omega(2) == 3 && par_omega(3) == 1
            damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
        elseif par_omega(2) == 1 && par_omega(3) == 0
            l = 4.3;
            damage.Sc = sqrt((1+1/3)*pi*G0*4.3/(8*E*horizon*0.66467));
        elseif par_omega(2) == 1 && par_omega(3) == 1
            l = 3;
            damage.Sc = sqrt(G0*(1/3+1)*pi^(3/2)/8/E*(l/horizon));
        else
            warning('Critical bond not defined.')
        end
        model.number = 1;
        % Damage dependent Sc
        if damage.DD
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        
    case "LBB"
        %% BB
        alfa = 0; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        mm = weightedVolume(par_omega); 
        c(1) = 6*E/mm;
        T = @models.forces.interactionForce_LLPSBB;
        model.linearity = true;
        model.stiffnessAnal = true; % true if an analytical stiffness matrix for such model is implemented
        model.dilatation = false;
        model.number = 2;
        
        
    case "LSJ"
        %% Lipton
        horizon = par_omega(1);
        alfa = 1;
        mm = weightedVolume(par_omega);
        %c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2; % Plane Strain
        %c(2) =(2*pi*horizon^3)^2/mm^2*E*1e6*(4*nu-1)/(2*(1+nu)*(1-2*nu))/2; %Plane Strain
        c(1) = 8*pi*horizon^3/mm*E/(1+nu)/2;
        c(2) = -2*(-1 + 3*nu)/(-1 + nu^2)*(pi*horizon^3/mm)^2*E;
        mu = E/(2*(1+nu));
        lambda = E*nu/(1+nu)/(1-2*nu);
        
        % -- Damage criteria
        if par_omega(2) == 3 && par_omega(3) == 1
            damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
        elseif par_omega(2) == 1 && par_omega(3) == 0
            damage.Sc = sqrt((1+nu)*pi*G0*4.3/(8*E*horizon*0.66467));
        else
            warning('Critical bond not defined.')
        end
        damage.thetaC = 3;
        
        T = @models.forces.interactionForce_Lipton;
        model.linearity = true;
        model.stiffnessAnal = false;
        model.dilatation = true;
        model.number = 3; 
        model.dilatHt = false;
        if damage.DD
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        
    case "LPS-T"
        %% LPS
        horizon = par_omega(1);
        alfa = 1;
        mm = weightedVolume(par_omega);
        kappa = E/3/(1-2*nu); mu = E/2/(1+nu);
        c(1) = kappa + mu/9*(nu+1)^2/(2*nu-1)^2;
        c(2) = 8*mu/mm;
        c(3) = nu;
        T = @models.forces.interactionForce_StateBased;
        model.linearity = false;
        model.stiffnessAnal = true;
        model.dilatation = true;
        % -- Sc
        if par_omega(2) == 3 && par_omega(3) == 1
            damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
        elseif par_omega(2) == 1 && par_omega(3) == 0
            l = 4.3;
            damage.Sc = sqrt((1+nu)*pi*G0*l/(8*E*horizon*0.66467));
        elseif par_omega(2) == 1 && par_omega(3) == 1
            l = 3;
            damage.Sc = sqrt(G0*(nu+1)*pi^(3/2)/8/E*(l/horizon));
        else
            warning('Critical bond not defined.')
        end
        % Damage dependent Sc
        if damage.DD
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        model.number = 4;
        
        
        case "PMB"
        %% PMB
        alfa = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(par_omega);        
        horizon = par_omega(1);
        c(1) = 6*E/mm;
        T = @models.forces.interactionForce_PMB;
        model.linearity = false;
        model.stiffnessAnal = true;
        model.dilatation = false;
        model.number = 5;
        
        % -- Sc
        if par_omega(2) == 3 && par_omega(3) == 1
            damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
        elseif par_omega(2) == 1 && par_omega(3) == 0
            l = 4.3;
            damage.Sc = sqrt((1+1/3)*pi*G0*l/(8*E*horizon*0.66467));
        elseif par_omega(2) == 1 && par_omega(3) == 1
            l = 3;
            damage.Sc = sqrt(G0*(1/3+1)*pi^(3/2)/8/E*(l/horizon));
        else
            warning('Critical bond not defined.')
        end
      
        % Damage dependent Sc
        if damage.DD
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        
        case "LSJ-T"
            horizon = par_omega(1);
            alfa = 1;
            mm = weightedVolume(par_omega);
            c(1) = 8*pi*horizon^3/mm*E/(1+nu)/2;
            c(2) = -2*(-1 + 3*nu)/(-1 + nu^2)*(pi*horizon^3/mm)^2*E;
            mu = E/(2*(1+nu));
            lambda = E*nu/(1+nu)/(1-2*nu);
            damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E)/horizon);
            damage.thetaC = 3;
            T = @models.forces.interactionForce_LSJT;
            model.linearity = true;
            model.stiffnessAnal = false;
            model.dilatation = true;
            model.number = 6; 
            model.dilatHt = true;
            % -- Sc
            if par_omega(2) == 3 && par_omega(3) == 1
                damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
            elseif par_omega(2) == 1 && par_omega(3) == 0
                l = 4.3;
                damage.Sc = sqrt((1+nu)*pi*G0*l/(8*E*horizon*0.66467));
            elseif par_omega(2) == 1 && par_omega(3) == 1
                l = 3;
                damage.Sc = sqrt(G0*(nu+1)*pi^(3/2)/8/E*(l/horizon));
            else
                warning('Critical bond not defined.')
            end
            % Damage dependent Sc
            if damage.DD
                damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
            else
                damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
            end
        case "PMB Concrete"
        %% PMB
        alfa = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(par_omega);
        ft = G0(1);
        fc = G0(2);
        horizon = par_omega(1);
        c(1) = 6*E/mm;
        T = @models.forces.interactionForce_PMBConcrete;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
        damage.St = ft/E;
        damage.Sc = -fc/E;
        model.number = 7;
        
        otherwise
        error("Chosen model is not implemented or it was mistyped");
end
end