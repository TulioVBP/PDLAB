function [ndof,idb,bc_set,bb,noFail] = boundaryCondition(x,stresses,m,h,A,option,tractionOpt,pc,damage)
%% INPUT parameter: 
% - x: the nodes position 
% - stresses: [sigma_x, sigma_y, tau_xy]
% - m: number of nodes inside delta
% - h: grid spacing
% - A: element area
% - option: boundary condition option (soon to be removed)
% - tractionOpt: where the stresses are applied
% - pc: prescribed constraints
% - damage: damage struct variable
%% OUTPUT parameter
% - ndof: Number of degree of freedoms
% - idb: collumn vector that index each of its row (corresponding to the
%        actual dof) to the matrix dof
%        EXAMPLE: idb = [1 3 4 6 7 8] means that the thirf row of the stiffness
%                matrix is related to the 2nd dof freedom of the system
% - bc_set: contains the set of constrained degrees of freedom on the first
%           collumn and its corresponding value on the second collumn; the
%           third collumn contains the corresponding dof velocity
% - bb: the actual body force vector
% - noFail: set of nodes for which we have no fail condition (mu = 1 always)

    %% DEFINE THE BOUNDARY CONDITIONS
    switch option
        case 1
        % Implement here your test
        % Example of BCs
        b = max(x(:,1));
        rightLay = find(x(:,1) == b);
        bottomLay = find(x(:,2) == min(x(:,2)));
        dof_Right = rightLay*2 - 1; % Constrain the x dof
        dof_Bot = bottomLay*2; % Constrain the y dof
        bc_set(:,1) = [dof_Right; dof_Bot]; % Each degree of freedom have a index number related
        bc_set(:,1) = sort(bc_set(:,1));
        bc_set(:,2:3) = zeros(length(bc_set(:,1)),2);
        B_given = [];
        case 0
        % Experiment itself
        a = max(x(:,2));
        bottomLay = find(x(:,2) == 0);
        topLay = find(x(:,2) == a);
        bc_set = [];
        B_given = [];
        case 2
            % 4 point bending test
            h = vecnorm(x(1,:) - x(2,:));
            b = max(x(:,1));
            bottomLay = find(x(:,2) == min(x(:,2)));
            topLay = find(x(:,2) == max(x(:,2)));
            bottom_constraint = bottomLay([2,length(bottomLay)-1]);
            top_constraint = topLay([2*floor(length(topLay)/5) 3*floor(length(topLay)/5)]);
            dof_bot = [bottom_constraint'*2 bottom_constraint(1)*2-1]';
            B_given = [top_constraint, zeros(length(top_constraint),1),-stresses(2)/h*ones(length(top_constraint),1)];
            bc_set(:,1) = [dof_bot]; % Each degree of freedom have a index number related
            bc_set(:,1) = sort(bc_set(:,1));
            bc_set(:,2:3) = zeros(length(bc_set(:,1)),2);
        case 3
            % Different prescribed case
            % - Displacement
            if ~isempty(pc.disp)
                dof_disp_constraint(:,1) = [2*(pc.disp(logical(pc.disp(:,2)),1))-1; 2*pc.disp(logical(pc.disp(:,3)),1)]; % [node boolean boolean value value]
                dof_disp_constraint(:,2) = [pc.disp(logical(pc.disp(:,2)),4); pc.disp(logical(pc.disp(:,3)),5)]; % Values
                bc_set(:,1) = dof_disp_constraint(:,1);
                bc_set(:,2) = dof_disp_constraint(:,2);
                bc_set(:,3) = zeros(size(bc_set(:,1)));
            end
            % - Velocity
            if ~isempty(pc.vel)
                dof_velocity_constraint(:,1) = [2*(pc.vel(logical(pc.vel(:,2)),1))-1; 2*pc.vel(logical(pc.vel(:,3)),1)]; % [boolean boolean value value]
                dof_velocity_constraint(:,2) = [pc.vel(logical(pc.vel(:,2)),4); pc.vel(logical(pc.vel(:,3)),5)]; % Values
                bc_set = [bc_set; dof_velocity_constraint(:,1),zeros(size(dof_velocity_constraint(:,1))), dof_velocity_constraint(:,2)];
            end
            
            if ~isempty(pc.vel) || ~isempty(pc.disp)
                [bc_set(:,1),II] = sort(bc_set(:,1)); % Sorting in ascending order
                bc_set(:,2:3) = bc_set(II,2:3); % Rearranging the displacement and velocity accordingly
            else
                bc_set = [];
            end
            % - Traction forces
            if ~isempty(pc.bodyForce)
                if length(size(pc.bodyForce)) == 2
                    % Constant body forces
                    B_given = pc.bodyForce(:,1:3); % [node, b_x, b_y]
                elseif length(size(pc.bodyForce)) == 3
                    B_given = pc.bodyForce(:,1:3,:); % [node, b_x, b_y, time]
                else
                    error("Wrong body force vector.")
                end
            else
                B_given = [];
            end
    end
    
    %% DEFINE THE IDB VECTOR
    ndof = 2*length(x) - length(bc_set);
    %in_ndof = 1:2*length(x);
    idb = zeros(2*length(x),1);
    id_dof = 1;
    id_const = ndof+1;
    for ii=1:2*length(x)
        if isempty(bc_set)
            check = 0;
        else
            check = sum(bc_set(:,1) == ii);
        end
         if check == 0
             % Free degree of freedom
             idb(ii) = id_dof;
             id_dof = id_dof + 1;
         else
             % Constrained degree of freedom
             idb(ii) = id_const;
             id_const = id_const + 1;
         end
    end
    % Ensure compatibility between idb for constrained nodes and bc_set -
    % TO BE DONE
    %% DEFINE THE BODY FORCE
    [b_old,noFail] = bodyForce(x,stresses, m, h, A,tractionOpt,B_given);
    if length(size(b_old)) == 2
        % Constant
        bb = zeros(2*length(x),1);
        for ii = 1:size(b_old,1)
            dofi = [idb(2*ii-1) idb(2*ii)];
            bb(dofi) = b_old(ii,:,:)'; % Assigning the body force values to the collumn vector
        end
    else
       % Non-Constant
        bb = zeros(2*length(x),size(b_old,3)); % 2Nxn matrix
        for ii = 1:size(b_old,1)
            dofi = [idb(2*ii-1) idb(2*ii)];
            bb(dofi,:) = permute(b_old(ii,:,:),[2,1,3]); % Assigning the body force values to the collumn vector
        end 
    end
   %% Updating no fail variable
   if ~isempty(bc_set)
    constraintNodes = sort(floor(bc_set(:,1)/2+2/3));
    noFail(constraintNodes) = 1;
   end
   if contains(tractionOpt,'off')
       noFail = zeros(size(noFail));
   end
    
   %% Plot the b.c.s
   figure 
   scatter(x(:,1),x(:,2),'b','filled','DisplayName','Free nodes')
   hold on
   if length(size(b_old)) == 2 && sum(sum(b_old ~= 0))
    scatter(x(b_old(:,1) ~= 0 | b_old(:,2)~=0,1),x(b_old(:,1) ~= 0 | b_old(:,2)~=0,2),'r','filled','DisplayName','Traction force nodes')
   elseif length(size(pc.bodyForce)) == 3
    scatter(x(b_old(:,1,end) ~= 0 | b_old(:,2,end)~=0,1),x(b_old(:,1,end) ~= 0 | b_old(:,2,end)~=0,2),'r','filled','DisplayName','Traction force nodes')
   end
   if ~isempty(bc_set)
    if ~isempty(bc_set(bc_set(:,3)== 0,:))
        scatter(x(floor(bc_set(bc_set(:,3)== 0,1)/2+2/3),1),x(floor(bc_set(bc_set(:,3)== 0,1)/2+2/3),2),'k','filled','DisplayName','Displacement const. nodes')
    end
    if ~isempty(bc_set(bc_set(:,3)~= 0,:))
        scatter(x(floor(bc_set(bc_set(:,3)~= 0,1)/2+2/3),1),x(floor(bc_set(bc_set(:,3)~= 0,1)/2+2/3),2),'g','filled','DisplayName','Velocity nodes')
    end
    %legend('Free nodes','Traction forces nodes','Displacement const. nodes','Velocity nodes')
   else
    %legend('Free nodes','Traction forces nodes')
   end
   if nargin > 8
      if ~isempty(damage.crackIn)
      x_crack = damage.crackIn(:,1);
      y_crack = damage.crackIn(:,2);
      line(x_crack,y_crack,'Color','red','LineWidth',4,'DisplayName','Initial crack')
      end
   end
   legend
   axis equal
   grid on
   xlabel('x (m)'); ylabel('y (m)')
   title('Boundary conditions')
   set(gca,'FontSize',13)
   pause(1)
end