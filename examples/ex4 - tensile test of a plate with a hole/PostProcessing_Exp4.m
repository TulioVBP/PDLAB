function PostProcessing_Exp4(x,u,n,idb,energy,phi,t,t_s)
close all
% Input: 
% - x = [x y]: position matrix
% - u: displacement matrix
% - n: evaluation time
% - phi: damage index
% - energy: struct variable with all the energies for the system
% - idb: come on, I know you know how this one works
% - t: simulation time vector
% - t_s: sampling time vector

%% Find n_sample
n_sample = find(t_s >= t(n),1);
if ~isempty(u)
    %% Making u a 3D matrix
    u = threeDModification(x,u,idb);
    %% Plot the displacement and strain map
    displacementPlot(x,u(:,:,n_sample),phi(:,n_sample));
    damagePlot(x,phi,n_sample);
    %strainPlot(x,u(:,:,n_sample)); % To be perfected
    %plotClassicalDerivative(x,u(:,:,n))
end

%% Plot the total energy
if exist('energy','var')~=0 
    if ~isempty(energy)
        energyPlot(x,energy,n_sample);
    end
end

end

function strainPlot(x,u)
    %% Function to plot the displacement
    v = u(:,1); % x component of the displacement
    w = u(:,2); % y component of the displacement field
    h = norm(x(1,:) - x(2,:)); % grid spacing
    x_lim = [min(x(:,1)) max(x(:,1))];
    y_lim = [min(x(:,2)) max(x(:,2))];
    bottom = find(x(:,2) < (y_lim(1) + 1e-12));
    top = find(x(:,2) > (y_lim(2) - 1e-12));
    left = find(x(:,1) < (x_lim(1) + 1e-12));
    right = find(x(:,1) > (x_lim(2) - 1e-12));
    N = 0; % Number of nodes in a row
    for ii = 1:length(x)
        if x(ii,2)~=x(1,2)
            break;
        end
        N = N+1;
    end
    out_layers = [top; bottom; left; right];
    for ii = 1:length(x)
%         if sum(ii == top) == 1 && sum(ii == right) == 1
%             % Top right corner
%             exx(ii) = (v(ii) - v(ii-1))/h; % dv/dx
%             eyy(ii) = (w(ii) - w(ii-N))/h; % dw/dy
%             exy(ii) = 1/2/h*(v(ii) - v(ii-N) + w(ii-1) - w(ii)); % 1/2*(dv/dy + dw/dx)
%         elseif sum(ii == top) == 1
%             % Top edge
%             exx(ii) = (v(ii+1) - v(ii))/h; % dv/dx
%             eyy(ii) = (w(ii) - w(ii-N))/h; % dw/dy
%             exy(ii) = 1/2/h*(v(ii) - v(ii-N) + w(ii+1) - w(ii)); % 1/2*(dv/dy + dw/dx)
%         elseif sum(ii == right) == 1
%             % Right edge
%             exx(ii) = (v(ii) - v(ii-1))/h; % dv/dx
%             eyy(ii) = (w(ii+N) - w(ii))/h; % dw/dy
%             exy(ii) = 1/2/h*(v(ii+N) - v(ii) + w(ii+1) - w(ii)); % 1/2*(dv/dy + dw/dx)
%         else
        if ~ismember(ii,out_layers)
            % Bulk nodes
            exx(ii) = (v(ii+1) - v(ii-1))/2/h; % dv/dx
            eyy(ii) = (w(ii+N) - w(ii-N))/2/h; % dw/dy
            exy(ii) = 1/4/h*(v(ii+N) - v(ii-N) + w(ii+1) - w(ii-1)); % 1/2*(dv/dy + dw/dx)
        else
            exx(ii) = 0; % dv/dx
            eyy(ii) = 0; % dw/dy
            exy(ii) = 0;
        end
    end
    %% Transforming into matrix
    h = norm(x(1,:) - x(2,:));
    b = max(x(:,1))-h/2 - min(x(:,1))+h/2;
    a = max(x(:,2))-h/2 - min(x(:,2))+h/2;
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
    EXX = zeros(size(X)); EYY = zeros(size(Y)); EXY = zeros(size(X));
    for jjj = 1:size(X,2)
        for iii = 1:size(Y,1)
            %ind = find(x(:,1) == X(1,jj) & x(:,2) == Y(ii,1)); not
            %reliable
            coord = [X(1,jjj) Y(iii,1)];
            ind = [];
            for ii = 1:length(x)
               if x(ii,1) < coord(1) + 1e-10 && x(ii,1) > coord(1) - 1e-6  && x(ii,2) < coord(2) + 1e-10 && x(ii,2) > coord(2) - 1e-10
                   ind = ii;
                   break;
               end
            end
            if ~isempty(ind)
                EXX(iii,jjj) = exx(ind);
                EYY(iii,jjj) = eyy(ind);
                EXY(iii,jjj) = exy(ind);
            else
                disp('Error: problem with the matching condition')
            end
        end        
    end
    %% Plot
    figure
    subplot(2,1,1)
    surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),EXX(2:end-1,2:end-1))
    xlabel x
    ylabel y
    zlabel \epsilon_x
    
    subplot(2,1,2)
    surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),EYY(2:end-1,2:end-1))
    xlabel x
    ylabel y
    zlabel \epsilon_y
end

function displacementPlot(x,u,phi)
    b = max(x(:,1)) - min(x(:,1));
    a = max(x(:,2)) - min(x(:,2));
    h = norm(x(1,:) - x(2,:));
    %% Plot displacement in the mesh
    scaleFactor = 1/5*abs(max(max(x(:,1))-min(x(:,1)),max(x(:,2))-min(x(:,2)))/max(max(u)));
        
    % Scatter
    figure
    scatter(x(:,1)+u(:,1)*scaleFactor,x(:,2)+u(:,2)*scaleFactor,[],phi,'filled')
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 0.4]);
    legend('Deformed configuration')
    xlabel x
    ylabel y
    set(gca,'FontSize',15)
    grid on
end

function damagePlot(x,phi,n)
    % Scatter
    figure
    for tt = 1:n
    scatter(x(:,1),x(:,2),[],phi(:,tt),'filled')
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 0.4]);
    legend('Damage configuration')
    xlabel x
    ylabel y
    set(gca,'FontSize',15)
    axis equal
    grid on
    title(strcat('t = ',num2str(n)))
    pause(0.001)
    end
end

function energyPlot(x,energy,n_final)
   %% Strain energy density
    b = max(x(:,1) - min(x(:,1)));
    a = max(x(:,2) - min(x(:,2)));
    h = norm(x(1,:) - x(2,:));
    % Transforming into a matrix array
    [X,Y] = meshgrid(min(x(:,1)):h:max(x(:,1))+1e-13, min(x(:,2)):h:max(x(:,2))+1e-13);
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
   %% Total energy evolution
   int = energy.W + energy.KE - energy.EW; % Total energy for each point
   t = 1:1:n_final; % Number of time steps
   figure
   plot(t,sum(energy.W(:,1:n_final)),'--','LineWidth',1.5)
   hold on
   plot(t,sum(energy.KE(:,1:n_final)),'-.','LineWidth',1.5)
   plot(t,sum(energy.EW(:,1:n_final)),'-','LineWidth',1.5)
   plot(t,sum(int(:,1:n_final)),'-','LineWidth',2.0)
   legend('Strain energy','Kinetic energy','External work','Total internal energy')
   xlabel('Time step n')
   ylabel('Energy (J/m)')
   set(gca,'FontSize',15)
   grid on

   
end


function u = threeDModification(x,u_2D,idb)
    if size(u_2D,1) == 2*length(x)
        u = zeros(length(x),2,size(u_2D,2));
        ii = (1:length(x))';
            dofi = [idb(2*ii-1) idb(2*ii)];
            for n = 1:size(u_2D,2)
                u_temp = u_2D(:,n);
                u(:,:,n) = u_temp(dofi);
            end
    else
        u = u_2D;
    end
    disp('Displacement vector converted');
end