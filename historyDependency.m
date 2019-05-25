% Function created to initialize history dependent variables
function history = historyDependency(x,maxNeigh,model)
% OUTPUT
% - history: 3D matrix with the variables

switch model.name
    case "PMB"
        history(:,:,1) = zeros(length(x),maxNeigh); % S_max
    case "Lipton Free Damage"
        history(:,:,1) = zeros(length(x),maxNeigh); % js integral
        history(:,:,2) = zeros(length(x),maxNeigh); % jtheta-x integral
        history(:,:,3) = zeros(length(x),maxNeigh); % jtheta-y integral
    otherwise
        history(:,:,1) = zeros(length(x),maxNeigh); % Don't have a meaning
end