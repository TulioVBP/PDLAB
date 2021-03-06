function w = influenceFunction(norma,parameters)%(norma,horizon,option)
%% INPUT:
% - norma: |xi|
% - parameters: [horizon option alpha]
%% OUTPUT
% - w: influence function
%global alpha
horizon = parameters(1); 
option = parameters(2); 
alfa = parameters(3);

switch option
    case 1 % Exponential
        l = horizon/3;
        w = exp(-(norma).^2/ (l^2))./norma.^alfa;
    case 2
        % P0
        w = 1./norma.^alfa; 
    case 3
        % P1
        w = (1-norma/horizon)./norma.^alfa; 
    case 4
        % P3
        w = (1 - 3*(norma/horizon).^2 + 2*(norma/horizon).^3)./norma.^alfa; 
    case 5
        % P5
        w = (1 - 10*(norma/horizon).^3 + 15*(norma/horizon).^4 - 6*(norma/horizon).^5)./norma.^alfa; 
    case 6
        % P7
        w = (1 - 35*(norma/horizon).^4 + 84*(norma/horizon).^5 - 70*(norma/horizon).^6 + 20*(norma/horizon).^7)./norma.^alfa;
    case 7
        % Singular
        w = (horizon./norma-1)./norma.^alfa;
end
% if w < 0 || norma > horizon
%     w = 0;
% end
w(w<0) = 0;
end
