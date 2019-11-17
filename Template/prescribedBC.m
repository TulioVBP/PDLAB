function pc = prescribedBC(x,stresses)
% OUTPUT
% - pc.disp/pc.vel = [node_ID boolean_X boolean_Y float_X float_Y]
%                node_ID: identify the node one is prescribing the displacement
%                boolean_X: true if x-component is prescribed
%                boolean_Y: true if y-component is prescribed
%                float_X: value of prescribed displacement/velocity in X
%                float_Y: value of prescribed displacement/velocity in Y
% - pc.bodyForce = [node_ID float_X float_Y];
% INPUT 
% The inputs are modificable as this function is customizable
% - x: nodes position matrix
% - stresses: can gather user-input tractions

    % DATA
    h = vecnorm(x(1,:)-x(2,:));
    if ~isempty(stresses)
        sigma_x = stresses(1);
        sigma_y = stresses(2);
        tau_xy = stresses(3);
    end
    % Modify here your boundary conditions
    if true
        % BC with prescribed body forces
        bottom_lay = find(x(:,2) < min(x(:,2))+1e-12);
        top_lay = find(x(:,2) > max(x(:,2))-1e-12);
        nodeSet1 = bottom_lay(2:5); l1 = length(nodeSet1);
        nodeSet2 =  bottom_lay(end-4:end-1); l2 = length(nodeSet2);
        pc.disp = [nodeSet1, ones(l1,1), ones(l1,1), zeros(l1,1), zeros(l1,1);
                   nodeSet2, zeros(l2,1), ones(l2,1), zeros(l2,1),zeros(l2,1)]; % [node boolean boolean value value]
        pc.vel = [];
        nodeSet3 = [top_lay(ceil(3/10*length(top_lay))-1:ceil(3/10*length(top_lay))+1); top_lay(ceil(7/10*length(top_lay))-1:ceil(7/10*length(top_lay))+1)];
        l3 = length(nodeSet3);
        pc.bodyForce = [nodeSet3,zeros(l3,1), -sigma_y/h*ones(l3,1)];
    else
        % BC with velocities instead - DYN ONLY (RIGID BODY)
        bottom_lay = find(x(:,2) < min(x(:,2))+1e-12);
        top_lay = find(x(:,2) > max(x(:,2))-1e-12);
        nodeSet1 = [bottom_lay(x(bottom_lay,1) > -0.45 - 1e-12 & x(bottom_lay,1) < -0.45+1e-12);...
                     bottom_lay(x(bottom_lay,1) > 0.45 - 1e-12 & x(bottom_lay,1) < 0.45+1e-12)];
        l1 = length(nodeSet1);
        pc.disp = [nodeSet1, zeros(l1,1), ones(l1,1), zeros(l1,1), zeros(l1,1)];
        nodeSet2 = [top_lay(x(top_lay,1) > -0.2 - 1e-12 & x(top_lay,1) < -0.2+1e-12);...
                     top_lay(x(top_lay,1) > 0.2 - 1e-12 & x(top_lay,1) < 0.2+1e-12)];
        vy = 0.1;
        pc.vel = [nodeSet2, zeros(l1,1), ones(l1,1), zeros(l1,1), -ones(l1,1)*vy];
        pc.bodyForce = [];
    end
end