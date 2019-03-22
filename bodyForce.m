function b = bodyForce(x,stresses,m,h,A)
%% INPUT:
% - x = [x y]: nodes' position vector
% - stresses = [sigma_x sigma_y tau_xy]: applied normal and shear stresses
% - m: number of spacings inside the horizon
% - h: grid spacing
% - A: element area vector
    x_lim = [min(x(:,1)) max(x(:,1))]; % Horizontal bounds
    y_lim = [min(x(:,2)) max(x(:,2))]; % Vertical bounds
    l = x_lim(2) - x_lim(1); % Length
    hh = y_lim(2) - y_lim(1); % Height
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
    %% Applying Neumann boundary conditions
    sigmax = stresses(1);
    sigmay = stresses(2);
    tauxy = stresses(3);
    b = zeros(size(x));
    % Top layer
    b(topLay,:) = b(topLay,:) + [tauxy,sigmay]*h/m./A(topLay);
    % Bottom layer
    b(bottomLay,:) = b(bottomLay,:) + [-tauxy,-sigmay]*h/m./A(bottomLay);
    % Left layer
    b(leftLay,:) = b(leftLay,:) + [-sigmax,-tauxy]*h/m./A(leftLay);
    % Right layer
    b(rightLay,:) = b(rightLay,:) + [sigmax,tauxy]*h/m./A(rightLay);
end