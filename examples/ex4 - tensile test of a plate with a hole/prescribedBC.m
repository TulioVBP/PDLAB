function pc = prescribedBC(x,stresses)
    % DATA
    h = vecnorm(x(1,:)-x(2,:));
    if ~isempty(stresses)
        sigma_x = stresses(1);
        sigma_y = stresses(2);
        tau_xy = stresses(3);
    end
    if false
    % Modify here your boundary conditions
        % BC with velocities instead - DYN ONLY (RIGID BODY)
        bottom_lay = find(x(:,2) < min(x(:,2))+1e-12);
        top_lay = find(x(:,2) > max(x(:,2))-1e-12);
        nodeSet1 = bottom_lay;%[bottom_lay(x(bottom_lay,1) > -0.45 - 1e-12 & x(bottom_lay,1) < -0.45+1e-12);...
                     %bottom_lay(x(bottom_lay,1) > 0.45 - 1e-12 & x(bottom_lay,1) < 0.45+1e-12)];
        l1 = length(nodeSet1);
        pc.disp = [nodeSet1, zeros(l1,1), ones(l1,1), zeros(l1,1), zeros(l1,1)];
        nodeSet2 = top_lay;%top_lay(x(top_lay,1) > -0.2 - 1e-12 & x(top_lay,1) < -0.2+1e-12);...
                     %top_lay(x(top_lay,1) > 0.2 - 1e-12 & x(top_lay,1) < 0.2+1e-12)];
        vy = 2;
        pc.vel = [nodeSet2, zeros(l1,1), ones(l1,1), zeros(l1,1), ones(l1,1)*vy];
        pc.bodyForce = [];
    else
        xc = [20-32.5,20-60;
              20-32.5,100-60]*10^-3;
        rc = [5,5]*10^-3;

        pin_up = find(vecnorm(x'-xc(2,:)')' < rc(2) + 3*h/4);
        pin_down = find(vecnorm(x'-xc(1,:)')' < rc(1) + 3*h/4);
        nodeSet1 = pin_down; l1 = length(nodeSet1);
        pc.disp = [nodeSet1, zeros(l1,1), ones(l1,1), zeros(l1,1), zeros(l1,1)];
        nodeSet2 = pin_up; l2 = length(nodeSet2);
        vy = 1;
        pc.vel = [nodeSet2, zeros(l2,1), ones(l2,1), zeros(l2,1), ones(l2,1)*vy];
        pc.bodyForce = [];
    end
end