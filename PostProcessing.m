function PostProcessing(x,u,n,idb,energy,phi)
% Input: 
% - x = [x y]: position matrix
% - family = family matrix
% - u: displacement matrix
% - n: evaluation time
% - phi: damage index
% - energy: struct variable with all the energies for the system
% - idb: come on, I know you know how this one works
%% Making u a 3D matrix
u = threeDModification(x,u,idb);
%% Plot the displacement and strain map
displacementPlot(x,u(:,:,n));
%strainPlot(x,u(:,(2*n-1):2*n)); % To be perfected
%% Plot the damage index
if exist('phi','var')~=0
    damagePlot(x,phi(:,n));
end
%% Plot the total energy
if exist('energy','var')~=0
    energyPlot(x,energy,n);
end
end



function strainPlot(x,u)
    %% Function to plot the displacement
    v = u(:,1); % x component of the displacement
    w = u(:,2); % y component of the displacement field
    h = norm(x(1,:) - x(2,:)); % grid spacing
    x_lim = [min(x(:,1)) max(x(:,1))];
    y_lim = [min(x(:,2)) max(x(:,2))];
    %bottom = find(x(:,2) < (y_lim(1) + 1e-12));
    top = find(x(:,2) > (y_lim(2) - 1e-12));
    %left = find(x(:,1) < (x_lim(1) + 1e-12));
    right = find(x(:,1) > (x_lim(2) - 1e-12));
    N = 0; % Number of nodes in a row
    for ii = 1:length(x)
        if x(ii,2)~=x(1,2)
            break;
        end
        N = N+1;
    end
    for ii = 1:length(x)
        if sum(ii == top) == 1 && sum(ii == right) == 1
            % Top right corner
            exx(ii) = (v(ii) - v(ii-1))/h; % dv/dx
            eyy(ii) = (w(ii) - w(ii-N))/h; % dw/dy
            exy(ii) = 1/2/h*(v(ii) - v(ii-N) + w(ii-1) - w(ii)); % 1/2*(dv/dy + dw/dx)
        elseif sum(ii == top) == 1
            % Top edge
            exx(ii) = (v(ii+1) - v(ii))/h; % dv/dx
            eyy(ii) = (w(ii) - w(ii-N))/h; % dw/dy
            exy(ii) = 1/2/h*(v(ii) - v(ii-N) + w(ii+1) - w(ii)); % 1/2*(dv/dy + dw/dx)
        elseif sum(ii == right) == 1
            % Right edge
            exx(ii) = (v(ii) - v(ii-1))/h; % dv/dx
            eyy(ii) = (w(ii+N) - w(ii))/h; % dw/dy
            exy(ii) = 1/2/h*(v(ii+N) - v(ii) + w(ii+1) - w(ii)); % 1/2*(dv/dy + dw/dx)
        else
            % Bulk nodes (or left/bottom nodes)
            exx(ii) = (v(ii+1) - v(ii))/h; % dv/dx
            eyy(ii) = (w(ii+N) - w(ii))/h; % dw/dy
            exy(ii) = 1/2/h*(v(ii+N) - v(ii) + w(ii+1) - w(ii)); % 1/2*(dv/dy + dw/dx)
        end        
    end
    % Plot
    figure
    plot3(x(:,1),x(:,2),exx,'.')
    xlabel x
    ylabel y
    zlabel exx
    
    figure
    plot3(x(:,1),x(:,2),eyy,'.')
    xlabel x
    ylabel y
    zlabel eyy
end

function displacementPlot(x,u)
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(0:h:b, 0:h:a);
    V = zeros(size(X)); W = zeros(size(Y));
    for jjj = 1:size(X,2)
        for iii = 1:size(Y,1)
            %ind = find(x(:,1) == X(1,jj) & x(:,2) == Y(ii,1)); not
            %reliable
            coord = [X(1,jjj) Y(iii,1)];
            for ii = 1:length(x)
               if x(ii,1) < coord(1) + 1e-10 && x(ii,1) > coord(1) - 1e-6  && x(ii,2) < coord(2) + 1e-10 && x(ii,2) > coord(2) - 1e-10
                   ind = ii;
                   break;
               end
            end
            if ~isempty(ind)
                V(iii,jjj) = u(ind,1);
                W(iii,jjj) = u(ind,2);
            else
                disp('Error: problem with the matching condition')
            end
        end        
    end
    figure
    surf(X,Y,V)
    xlabel x
    ylabel y
    zlabel ux
    set(gca,'FontSize',15)
    
    figure
    surf(X,Y,W)
    xlabel x
    ylabel y
    zlabel uy
    set(gca,'FontSize',15)
    %% Plot displacement in the mesh
    scaleFactor = 1e3;
    % Mesh
    figure
    mesh(X,Y,zeros(size(X)),'LineWidth',1)
    hold on
    colormap winter
    mesh(X+scaleFactor*V,Y+scaleFactor*W,20*ones(size(X)),'LineWidth',1)
    hidden off
    xlabel x
    ylabel y
    set(gca,'FontSize',15)
    legend('Reference','Deformed')
%     % Scatter
%     figure
%     scatter(x(:,1),x(:,2),'b','filled')
%     hold on
%     scatter(x(:,1)+u(:,1)*scaleFactor,x(:,2)+u(:,2)*scaleFactor,'r','filled')
%     legend('Reference configuration','Deformed configuration')
%     xlabel x
%     ylabel y
%     set(gca,'FontSize',15)
%     grid on
end

function damagePlot(x,phi)
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(0:h:b, 0:h:a);
    PHI = zeros(size(X));
    for jjj = 1:size(X,2)
        for iii = 1:size(Y,1)
            %ind = find(x(:,1) == X(1,jj) & x(:,2) == Y(ii,1)); not
            %reliable
            coord = [X(1,jjj) Y(iii,1)];
            for ii = 1:length(x)
               if x(ii,1) < coord(1) + 1e-10 && x(ii,1) > coord(1) - 1e-6  && x(ii,2) < coord(2) + 1e-10 && x(ii,2) > coord(2) - 1e-10
                   ind = ii;
                   break;
               end
            end
            if ~isempty(ind)
                PHI(iii,jjj) = phi(ind);
            else
                disp('Error: problem with the matching condition')
            end
        end        
    end
    figure
    pcolor(X,Y,PHI)
    xlabel x
    ylabel y
    title('Damage index')
    set(gca,'FontSize',15)
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 1]);
end

function energyPlot(x,energy,n_final)
   %% Total energy evolution
   int = energy.W + energy.KE - energy.EW; % Total energy for each point
   t = 1:1:n_final; % Number of time steps
   figure
   plot(t,sum(energy.W(:,1:n_final)),'--','LineWidth',1.5)
   hold on
   plot(t,sum(energy.KE(:,1:n_final)),'-.','LineWidth',1.5)
   plot(t,sum(energy.EW(:,1:n_final)),'-d','LineWidth',1.5)
   plot(t,sum(int(:,1:n_final)),'-','LineWidth',1.5)
   legend('Strain energy','Kinetic energy','External work','Total internal energy')
   xlabel('Time step n')
   ylabel('Energy (J/m)')
   set(gca,'FontSize',15)
   grid on
   %% Strain energy density
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(0:h:b, 0:h:a);
    WW = zeros(size(X));
    for jjj = 1:size(X,2)
        for iii = 1:size(Y,1)
            %ind = find(x(:,1) == X(1,jj) & x(:,2) == Y(ii,1)); not
            %reliable
            coord = [X(1,jjj) Y(iii,1)];
            for ii = 1:length(x)
               if x(ii,1) < coord(1) + 1e-10 && x(ii,1) > coord(1) - 1e-6  && x(ii,2) < coord(2) + 1e-10 && x(ii,2) > coord(2) - 1e-10
                   ind = ii;
                   break;
               end
            end
            if ~isempty(ind)
                WW(iii,jjj) = energy.W(ind,n_final);
            else
                disp('Error: problem with the matching condition')
            end
        end        
    end
    figure
    surf(X,Y,WW)
    xlabel x
    ylabel y
    title('Strain Energy Density')
    set(gca,'FontSize',15)
   
end

function u = threeDModification(x,u_2D,idb)
    if size(u_2D,1) == 2*length(x)
        u = zeros(length(x),2,size(u_2D,2));
        for ii = 1:length(x)
            dofi = [idb(2*ii-1) idb(2*ii)];
            for n = 1:size(u_2D,2)
                u(ii,:,n) = u_2D(dofi,n)';
            end
        end
    else
        u = u_2D;
    end
    disp('Displacement vector converted');
end