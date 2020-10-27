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
        modelo = models.modelLSJ(E,nu,par_omega,damage,G0,dt);
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