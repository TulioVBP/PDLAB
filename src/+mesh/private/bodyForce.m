function [b,noFailZone] = bodyForce(x,stresses,m,h,A,tractionOpt,B_given)
% NEEDED UPDATE: Make it more general so that one can apply non-homogeneous
%                traction forces
%% INPUT:
% - x = [x y]: nodes' position vector
% - stresses = [sigma_x sigma_y tau_xy]: applied normal and shear stresses
% - m: number of spacings inside the horizon
% - h: grid spacing
% - A: element area vector or scalar
%% OUTPUT:
% b - body in N*2D condition
% noFailZone - set of nodes in the boundary layer of the non-null applied stresses
if length(A) == 1
    A = A*ones(length(x),1);
end
    x_lim = [min(x(:,1)) max(x(:,1))]; % Horizontal bounds
    y_lim = [min(x(:,2)) max(x(:,2))]; % Vertical bounds
    l = x_lim(2) - x_lim(1); % Length
    hh = y_lim(2) - y_lim(1); % Height
    % Correcting the limits
    x_lim = [min(x(:,1))-h/2 max(x(:,1))+h/2]; % Horizontal bounds
    y_lim = [min(x(:,2))-h/2 max(x(:,2))+h/2]; % Vertical bounds
    %% Defining the edges indexes
%     bottom = find(x(:,2) == y_lim(1)); % Bottom edge nodes index
%     top = find(x(:,2) == y_lim(2)); % Top edge nodes index
%     left = find(x(:,1) == x_lim(1)); % Left edge nodes index
%     right = find(x(:,1) == x_lim(2)); % Right edge nodes index
    %% Defining the boundary layers indexes
    bottomLay = find(x(:,2) < (y_lim(1) + m*h + 1e-12)); % Bottom layer nodes index
    topLay = find(x(:,2) > (y_lim(2) - m*h - 1e-12)); % Top layer nodes index
    leftLay = find(x(:,1) < (x_lim(1) + m*h + 1e-12)); % Left layer nodes index
    rightLay = find(x(:,1) > (x_lim(2) - m*h - 1e-12)); % Right layer nodes index
    % Layers for no fail zone
    bottomLayer = find(x(:,2) < (y_lim(1) + 3*m*h + 1e-12)); % Bottom layer nodes index
    topLayer = find(x(:,2) > (y_lim(2) - 3*m*h - 1e-12)); % Top layer nodes index
    %leftLayer = find(x(:,1) < (x_lim(1) + 3*m*h + 1e-12)); % Left layer nodes index
    %rightLayer = find(x(:,1) > (x_lim(2) - 3*m*h - 1e-12)); % Right layer nodes index
    %% Applying constant Neumann boundary conditions
    if ~isempty(stresses)
        sigmax = stresses(1);
        sigmay = stresses(2);
        tauxy = stresses(3);
    end
    b = zeros(size(x));
    
    if contains(tractionOpt,'-') && ~isempty(B_given)
       if length(size(B_given)) == 2
           % Constant
            b(B_given(:,1),:) = b(B_given(:,1),:) + B_given(:,2:3);
       else
           b = zeros([size(x) size(B_given,3)]);
           % Non-constant
           b(B_given(:,1),:,:) = b(B_given(:,1),:,:) + B_given(:,2:3,:);
       end
    else
        if contains(tractionOpt,'t')
            % Top layer
            b(topLay,:) = b(topLay,:) + [tauxy,sigmay]*h/m./A(topLay);
        end
        if contains(tractionOpt,'b')
            % Bottom layer
            b(bottomLay,:) = b(bottomLay,:) + [-tauxy,-sigmay]*h/m./A(bottomLay);
        end
        if contains(tractionOpt,'l')
            % Left layer
            b(leftLay,:) = b(leftLay,:) + [-sigmax,-tauxy]*h/m./A(leftLay);
        end
        if contains(tractionOpt,'r')
            % Right layer
            b(rightLay,:) = b(rightLay,:) + [sigmax,tauxy]*h/m./A(rightLay);
        end
    end
    %% Defining no fail zone
    noFailZone = zeros(length(x),1);
    if length(size(b)) == 3
        noFailZone(b(:,1,end)~=0|b(:,2,end)~=0) = 1;
    else
        noFailZone(b(:,1)~=0|b(:,2)~=0) = 1;
    end
    %noFailZone(topLayer) = true;
    %noFailZone(bottomLayer) = true;
    %noFailZone = find(b(:,1)~= 0 | b(:,2)~= 0); % Criterion 1: if a traction force is applied to the node, it is no fail zone
end