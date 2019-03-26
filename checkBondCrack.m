% Function to check if the bond is crossing a specific line
function [P,check] = checkBondCrack(x1, x2, x3, x4)
%% INPUT
% - x1 = [x y]: is the initial node
% - x2 = [x y]: is the final node
% - x3 = [x,y]: is the initial point of the crack
% - x4 = [x,y]: is the final point of the crack
%% OUTPUT
% - P =[Px Py]: is the intersection point
% - check = true or false: is true if the intersection point belong to both
%                          the bond and to the crack segment

den = det([det([x1(1) 1;x2(1) 1]) det([x1(2) 1;x2(2) 1]);det([x3(1) 1;x4(1) 1]) det([x3(2) 1;x4(2) 1])]);
%% X-component
det11 = det([x1(1) x1(2); x2(1) x2(2)]);
det12 = det([x1(1) 1; x2(1) 1]);
det21 = det([x3(1) x3(2); x4(1) x4(2)]);
det22 = det([x3(1) 1; x4(1) 1]);
Px = det([det11 det12; det21 det22])/den ;

%% Y-component
%det11 = det([x1(1) x1(2); x2(1) x2(2)]);
det12 = det([x1(2) 1; x2(2) 1]);
%det21 = det([x3(1) x3(2); x4(1) x4(2)]);
det22 = det([x3(2) 1; x4(2) 1]);
Py = det([det11 det12; det21 det22])/den ;

P = [Px Py];
if norm(P) ~= Inf
    %% Check if P belongs to both the crack segment and the bond
    if P(1) > min(x3(1),x4(1))-1e-12 && P(1) < max(x3(1),x4(1))+1e-12 && P(2) > min(x3(2),x4(2))-1e-12 && P(2) < max(x3(2),x4(2))+1e-12
        % Intersection belong to the crack segment
        if P(1) > min(x1(1),x2(1))-1e-12 && P(1) < max(x1(1),x2(1))+1e-12 && P(2) > min(x1(2),x2(2))-1e-12 && P(2) < max(x1(2),x2(2))+1e-12
            % Intersection point belong to the bond
            check = true;
        else
            % Intersection doesn't belong to the bond
            check = false;
        end
    else
        % Intersection doesn't belong to the crack segment
        check = false;
    end
else
    %% Infinite norm means parallel lines
    check = false;
end



end