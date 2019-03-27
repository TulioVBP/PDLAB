function PostProcessing(x,u,n,phi)
% Input: 
% - x = [x y]: position matrix
% - u: displacement matrix
% - n: evaluation time
% - phi: damage index

%% Plot the displacement and strain map
displacementPlot(x,u(:,(2*n-1):2*n));
%strainPlot(x,u(:,(2*n-1):2*n)); % To be perfected
%% Plot the damage index
damagePlot(x,phi(:,n)); 
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
    global a b
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
    
    figure
    surf(X,Y,W)
    xlabel x
    ylabel y
    zlabel uy
end

function damagePlot(x,phi)
    global a b
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
    c = jet(1000);
    colormap(c);
    colorbar
    caxis([0 1]);
end