function [modelo,damage] = modelParameters(model,par_omega,damage,E,nu,G0,dt)
% INPUTS
% - model: the constitutive model you are using. Implemented models are:
% "DTT", "PMB", "LBB", "LSJ", "LSJ-T","LPS-T"
% - par_omega = [horizon omega alfa]
% - damage: struct variable with damage parameters
% - E: Young modulus
% - nu: Poisson's ratio
% - G0: Energy release rate
% - dt: for LSJT model only
% OUTPUTS
% - modelo: the object with the constitutive model you are considering
% - damage: updated struct variable (soon to be deprecated)

switch model.name
    case "DTT"
    %% PMB DTT
        modelo = models.modelDTT(E,par_omega,damage,G0);
        
    case "LBB"
    %% Linear Bond-Based
        modelo = models.modelLBB(E,par_omega,damage,G0);
        
    case "LPS-T"
    %% LPS
        modelo = models.modelLPST(E,nu,par_omega,damage,G0);        

    case "PMB"
    %% PMB
        modelo = models.modelPMB(E,par_omega,damage,G0);

    case "LSJ-T"
        modelo = models.modelLSJT(E,nu,par_omega,damage,G0,dt);
        
    case "LSJ-T2"
        modelo = models.modelLSJT2(E,nu,par_omega,damage,G0,dt);

    case "LSJ"
            %% Lipton
            horizon = par_omega(1);
            modelo = models.modelLPST(E,nu,par_omega);

            % -- Damage criteria
            if par_omega(2) == 3 && par_omega(3) == 1
                damage.Sc = sqrt(5*pi*G0/9/(E)/horizon);
            elseif par_omega(2) == 1 && par_omega(3) == 0
                damage.Sc = sqrt((1+nu)*pi*G0*4.3/(8*E*horizon*0.66467));
            else
                warning('Critical bond not defined.')
            end
            damage.thetaC = 3;

    %         T = @models.forces.interactionForce_Lipton;
    %         model.linearity = true;
    %         model.stiffnessAnal = false;
    %         model.dilatation = true;
    %         model.number = 3; 
    %         model.dilatHt = false;
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