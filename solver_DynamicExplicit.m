% Explicit time - dynamic solver

function [u_n,phi] = solver_DynamicExplicit(x,maxNeigh,t,idb,b,familyMat,partialAreas,T,crackIn)
    global rho
    dt = abs(t(2) - t(1));
    % ########### INITIALIZE SIMULATION MATRICES ##############
        % M = rho * eye(length(x)); % Inertia matrix
        % Minv = 1/rho * eye(length(x)); % Inverse of the inertial matrix
        Minv = 1/rho; % Diagonal and with the same value: scalar
        S_max = zeros(length(x),maxNeigh);
        phi = zeros(length(x),length(t)); % No initial damage
        % Initial condition
        u_n = zeros(2*length(x),length(t)); % 2*N x n  2D matrix
        v_n = zeros(2*length(x),3); % Velocity Verlet matrix [i i+1/2]: by initializing it to zero, the rigid motion is eliminated.ffz
        bn = b;%bodyForce(x,stresses,m,h,A);
        for n = 1:length(t)-1
            % TIME LOOP
            fn = zeros(2*length(x),1);
            for ii = 1:length(x)
                dofi = [idb(2*ii-1) idb(2*ii)];
                % Loop on the nodes
                family = familyMat(ii,familyMat(ii,:)~=0);
               for jj = family
                   dofj = [idb(2*jj-1) idb(2*jj)];
                   % Loop on their neighbourhood
                   neig_index = find(familyMat(ii,:) == jj);
                   [fij,S_max(ii,neig_index)] = T(x(ii,:),x(jj,:),u_n(dofi,n)',u_n(dofj,n)',S_max(ii,neig_index),crackIn);
                   Vj = partialAreas(ii,neig_index);
                   fn(dofi) = fn(dofi) + (fij')*Vj;   
               end
               %phi(ii,n) = damageIndex(x,u_n(:,(2*n-1):(2*n)),family(ii,:),partialAreas(ii,:),ii,crackIn); % Damage index
               phi(ii,n) = damageIndex(x,u_n(:,n),familyMat(ii,:),partialAreas(ii,:),ii,crackIn,idb); % Damage index
            end
            % ############ VELOCITY VERLET ALGORITHM ###############
            %a_n = dt/2*M\(fn + bn); % Acceleration
            % Step 1 - Midway velocity
            v_n(:,2) = v_n(:,1) + dt/2*Minv*(fn + bn); % V(n+1/2)
            %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
            u_n(:,n+1) = u_n(:,n) + dt*v_n(:,2); % u(n+1)
            v_n(:,1) = v_n(:,2) + dt/2*Minv*(fn + bn); % V(n+1) is stored in the next V(n)
            % ############ COUNTING THE PROCESSING TIME #############
            disp("Percentage of the process: " + num2str(n/(length(t)-1)*100) + "%")
            %pause
        end
end
% 
%         % ########### INITIALIZE SIMULATION MATRICES ##############
%         % M = rho * eye(length(x)); % Inertia matrix
%         % Minv = 1/rho * eye(length(x)); % Inverse of the inertial matrix
%         Minv = 1/rho; % Diagonal and with the same value: scalar
%         S_max = zeros(length(x),maxNeigh);
%         phi = zeros(length(x),length(t)); % No initial damage
%         % Initial condition
%         u_n = zeros(length(x),2,length(t)); % N x 2 x n  3D matrix
%         %u = zeros(size(x));
%         v_n = zeros(length(x),2*3); % Velocity Verlet matrix [i i+1/2 i+1]: by initializing it to zero, the rigid motion is eliminated.ffz
%         bn = bodyForce(x,stresses,m,h,A);
%         for n = 1:length(t)-1
%             %% TIME LOOP
%             fn = zeros(size(x));
%             for ii = 1:length(x)
%                 % Loop on the nodes
%                for neig_index=1:size(family,2)
%                    % Loop on their neighbourhood
%                    jj = family(ii,neig_index); 
%                    if jj == 0 % Check if the jj neighbor is a valid node
%                        break
%                    elseif ii~=jj % Check if jj is a different node (shouldn't be necessary as family already accounts for it)
%                        % Bond doesn't cross the crack initial segment
%                        %[fij,S_max(ii,neig_index)] = interactionForce(x(ii,:),x(jj,:),u_n(ii,(2*n-1):(2*n)),u_n(jj,(2*n-1):(2*n)),S_max(ii,neig_index),crackIn);
%                        switch model
%                            case "PMB"
%                            [fij,S_max(ii,neig_index)] = interactionForce_PMB(x(ii,:),x(jj,:),u_n(ii,:,n),u_n(jj,:,n),S_max(ii,neig_index),crackIn);
%                            case "Linearized LPS bond-based"
%                            [fij,S_max(ii,neig_index)] = interactionForce_LLPSBB(x(ii,:),x(jj,:),u_n(ii,:,n),u_n(jj,:,n),S_max(ii,neig_index),crackIn);
%                            otherwise
%                            disp("Chosen model is not implemented or it was mistyped");
%                            quit
%                        end
%                            Vj = partialAreas(ii,neig_index);
%                            fn(ii,:) = fn(ii,:) + fij*Vj;
%                    end
%                end
%                %phi(ii,n) = damageIndex(x,u_n(:,(2*n-1):(2*n)),family(ii,:),partialAreas(ii,:),ii,crackIn); % Damage index
%                phi(ii,n) = damageIndex(x,u_n(:,:,n),family(ii,:),partialAreas(ii,:),ii,crackIn); % Damage index
%             end
%             % ############ VELOCITY VERLET ALGORITHM ###############
%             %a_n = dt/2*M\(fn + bn); % Acceleration
%             % Step 1 - Midway velocity
%             v_n(:,3:4) = v_n(:,1:2) + dt/2*Minv*(fn + bn); % V(n+1/2)
%             %u_n(:,(2*(n+1)-1):2*(n+1)) = u_n(:,(2*n-1):2*n) + dt*v_n(:,3:4); % u(n+1)
%             u_n(:,:,n+1) = u_n(:,:,n) + dt*v_n(:,3:4); % u(n+1)
%             v_n(:,5:6) = v_n(:,3:4) + dt/2*Minv*(fn + bn); % V(n+1) - unnecessary
%             v_n(:,1:2) = v_n(:,5:6); % For the next iterative loop
%             % ############ COUNTING THE PROCESSING TIME #############
%             disp("Stress: "+ int2str(s_index) + " - Mesh: " +int2str(m_index) + " - Percentage of the process: " + num2str(n/(length(t)-1)*100) + "%")
%             %pause