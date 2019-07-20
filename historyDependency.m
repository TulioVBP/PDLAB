% Function created to initialize history dependent variables
function history = historyDependency(x,maxNeigh,model)
% OUTPUT
% - history: 3D matrix with the variables

switch model.number
    case 1%"PMB"
            history.S = zeros(length(x),maxNeigh); % S_max
    case 5%"PMB"
            history.S = zeros(length(x),maxNeigh); % S_max
    case 3%"Lipton Free Damage"
            history.S = zeros(length(x),maxNeigh); % js integral
            history.theta = zeros(length(x),1); % jtheta-x integral
        %history(:,:,3) = zeros(length(x),maxNeigh); % jtheta-y integral
    case 4 %"LPS 2D"
            history.S = zeros(length(x),maxNeigh); % js integral
            history.theta = zeros(length(x),1); % jtheta-x integral
    otherwise
        history.S = zeros(length(x),maxNeigh); % Don't have a meaning
end