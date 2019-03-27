function mu = damageFactor(x,notch,x_i,x_j)
%% Input
% x: stretch
% crack: coordinates of the beginning and end points
% x_i: position of node i
% x_j: position of node j
%% Output
% mu: damage factor
%% CODE
    global S0 S1
    crackSegments = size(notch,1); % At least 2
    for ii = 1:crackSegments-1
        [~,check(ii)] = checkBondCrack(x_i,x_j,notch(ii,:),notch(ii+1,:));
    end
    if ~sum(check)
        % Bond doesn't interceptate the initial crack (notch)
        if x < S0(2) && x > -1
            mu = 1;
        elseif x > S0(2) && x < S1(2)
            mu = (S1(2) - x)/(S1(2) - S0(2));
        else
            mu = 0;
        end
    else
        % Bond interceptate the initial crack (notch)
      mu = 0;
    end
end