function pc = prescribedBC(x,stresses,d)
    if nargin < 3
        d = 4;
    end
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
        bottom_lay = find(x(:,2) < (min(x(:,2))-h/2 + 1*h +1e-12));
        top_lay = find(x(:,2) > (max(x(:,2))+h/2- 1*h - 1e-12));
        nodeSet1 = bottom_lay; l1 = length(nodeSet1);
        nodeSet2 =  top_lay; l2 = length(nodeSet2);
        pc.disp = []; % [node boolean boolean value value]
        pc.vel = [];
        % Sigma function
        sigma_fun = @(t,t0) (t<t0).*sigma_y.*t/t0 + (t>=t0).*sigma_y;
        t0 = 50;
        for n = 1:1301
            pc.bodyForce(:,:,n) = [nodeSet1,zeros(l1,1), -sigma_fun(65/1301*n,t0)/h/1*ones(l1,1);
                            nodeSet2,zeros(l2,1), sigma_fun(65/1301*n,t0)/h/1*ones(l1,1)];
        end
    else
        % BC with prescribed body forces constants
        bottom_lay = find(x(:,2) < (min(x(:,2))-h/2+ d*h+1e-12));
        top_lay = find(x(:,2) > (max(x(:,2))+h/2- d*h - 1e-12));
        nodeSet1 = bottom_lay; l1 = length(nodeSet1);
        nodeSet2 =  top_lay; l2 = length(nodeSet2);
        pc.disp = []; % [node boolean boolean value value]
        pc.vel = [];
        pc.bodyForce = [nodeSet1,zeros(l1,1), -sigma_y/h/d*ones(l1,1);
                        nodeSet2,zeros(l2,1), sigma_y/h/d*ones(l2,1)];
    end
end