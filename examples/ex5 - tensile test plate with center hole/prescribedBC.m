function pc = prescribedBC(x,horizon)
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

    % DATA
    h = vecnorm(x(1,:)-x(2,:));
    % Modify here your boundary conditions
    % BC with prescribed displacement
    left_lay = find(x(:,1) < min(x(:,1)) + horizon + 1e-12);
    right_lay = find(x(:,1) > max(x(:,1)) - horizon - 1e-12);
    mid_line = left_lay(x(left_lay,2) > - h/2 + 1e-12 & x(left_lay,2) < h/2 -1e-12);
    if isempty(mid_line)
        warning('Rigid displacement BC not prescribed')
    end
    %left_lay = left_lay(x(left_lay,2) < - h/2 + 1e-12 | x(left_lay,2) > h/2 -1e-12);
       
    nodeSet1 = left_lay; l1 = length(nodeSet1);
    nodeSet2 =  right_lay; l2 = length(nodeSet2);
    nodeSet3 = mid_line; l3 = length(nodeSet3);
    
    pc.disp = [nodeSet1, ones(l1,1), zeros(l1,1), zeros(l1,1), zeros(l1,1);
               nodeSet2, ones(l2,1), zeros(l2,1), ones(l2,1),zeros(l2,1);
               nodeSet3, ones(l3,1), ones(l3,1), zeros(l3,1), zeros(l3,1)]; % [node boolean boolean value value]
    pc.vel = [];
    pc.bodyForce = [];
end