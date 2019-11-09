function [model,c,T,damage] = modelParameters(model,par_omega,damage,E,nu,G0)
    switch model.name
    case "PMB DTT"
        %% PMB DTT
        alfa = 1; % Because for the PMB we always have to modulate the influence function by 1/|\xi|
        mm = weightedVolume(par_omega);
        horizon = par_omega(1);
        c(1) = 6*E*1e6/mm;
        T = @interactionForce_PMBDTT;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
        damage.Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon);
        model.number = 1;
        % Damage dependent Sc
        if false
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        
    case "Linearized LPS bond-based"
        %% BB
        alfa = 0; % If alpha = 1, I'm basically reducing the linearized bond based model to the linearized PMB model
        mm = weightedVolume(par_omega); 
        c(1) = 6*E*1e6/mm;
        T = @interactionForce_LLPSBB;
        model.linearity = true;
        model.stiffnessAnal = true; % true if an analytical stiffness matrix for such model is implemented
        model.dilatation = false;
        model.number = 2;
        
    case "Lipton Free Damage"
        %% Lipton
        horizon = par_omega(1);
        alfa = 1;
        mm = weightedVolume(par_omega);
        %c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2; % Plane Strain
        %c(2) =(2*pi*horizon^3)^2/mm^2*E*1e6*(4*nu-1)/(2*(1+nu)*(1-2*nu))/2; %Plane Strain
        c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2;
        c(2) = -2*(-1 + 3*nu)/(-1 + nu^2)*(pi*horizon^3/mm)^2*E*1e6;
        mu = E*1e6/(2*(1+nu));
        lambda = E*1e6*nu/(1+nu)/(1-2*nu);
        damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E*1e6)/horizon);
        damage.thetaC = 3;
        T = @interactionForce_Lipton;
        model.linearity = true;
        model.stiffnessAnal = false;
        model.dilatation = true;
        model.number = 3; 
        model.dilatHt = false;
        
    case "LPS 2D"
        %% LPS
        horizon = par_omega(1);
        alfa = 1;
        mm = weightedVolume(par_omega);
        kappa = E*1e6/3/(1-2*nu); mu = E*1e6/2/(1+nu);
        c(1) = kappa + mu/9*(nu+1)^2/(2*nu-1)^2;
        c(2) = 8*mu/mm;
        c(3) = nu;
        T = @interactionForce_StateBased;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = true;
        damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E*1e6)/horizon);
        % Damage dependent Sc
        if false
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
        c(1) = 6*E*1e6/mm;
        T = @interactionForce_PMB;
        model.linearity = false;
        model.stiffnessAnal = false;
        model.dilatation = false;
        damage.Sc = sqrt(5*pi*G0/9/(E*1e6)/horizon);
        model.number = 5;
        % Damage dependent Sc
        if true
            damage.alfa = 0.2; damage.beta = 0.2; damage.gamma = 1.4;
        else
            damage.alfa = 0; damage.beta = 0; damage.gamma = 1; % No dependency
        end
        
        case "LSJ-T"
            horizon = par_omega(1);
            alfa = 1;
            mm = weightedVolume(par_omega);
            %c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2; % Plane Strain
            %c(2) =(2*pi*horizon^3)^2/mm^2*E*1e6*(4*nu-1)/(2*(1+nu)*(1-2*nu))/2; %Plane Strain
            c(1) = 8*pi*horizon^3/mm*E*1e6/(1+nu)/2;
            c(2) = -2*(-1 + 3*nu)/(-1 + nu^2)*(pi*horizon^3/mm)^2*E*1e6;
            mu = E*1e6/(2*(1+nu));
            lambda = E*1e6*nu/(1+nu)/(1-2*nu);
            damage.Sc = sqrt(5*(1+nu)*pi*G0/12/(E*1e6)/horizon);
            damage.thetaC = 3;
            T = @interactionForce_LSJT;
            model.linearity = true;
            model.stiffnessAnal = false;
            model.dilatation = true;
            model.number = 6; 
            model.dilatHt = true;
            
        otherwise
        error("Chosen model is not implemented or it was mistyped");
end
end