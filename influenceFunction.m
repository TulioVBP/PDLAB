function w = influenceFunction(norma,horizon,option)
%% INPUT:
% - norma: |xi|
% - horizon: horizon
% - option: chosen influence function
%% OUTPUT
% - w: influence function
global alpha

switch option
    case 1
        l = horizon/4.3;
        w = exp(-(norma)^2/ (l^2));
    case 2
        % P0
        w = 1/norma^alpha; 
    case 3
        % P1
        w = (1-norma/horizon)/norma^alpha; 
    case 4
        % P3
        w = (1 - 3*(norma/horizon)^2 + 2*(norma/horizon)^3)/norma^alpha; 
    case 5
        % P5
        w = (1 - 10*(norma/horizon)^3 + 15*(norma/horizon)^4 - 6*(norma/horizon)^5)/norma^alpha; 
    case 6
        % P7
        w = (1 - 35*(norma/horizon)^4 + 84*(norma/horizon)^5 - 70*(norma/horizon)^6 + 20*(norma/horizon)^7)/norma^alpha;
end
% if w < 0
%     w = 0;
% end
end
