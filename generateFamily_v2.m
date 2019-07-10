function [family,partialAreas,maxNeigh,surfaceCorrection] = generateFamily_v2(x,horizon,d,option,test,PA_alg,surf_Alg)
% PA-HHB IMPLEMENTATION
% Volume correction implementation
% INPUT : 
%---------
% x - Position of the nodes
% horizon - Peridynamic horizon
% d - number of intervals taken by horizon
% test - boolean variable that takes true if we are performing a test on
% the script
% PA_alg - Partial area algorithm
% ------
% OUTPUT:
%---------
% family - Each line 'i' of the family is related to the node of position
%          x(i,:), until it becomes zero
% familyInfo - It informs the A_k^(i) (partial area) whether the node k in the neighborhood i is on the
%              border of the neighborhood or not;
% maxNeigh - the maximum number of neighbors around a node (scalar)
% ------
%% Entering in the directory
if test
    dirname = 'FAMILY FILES TEST';
else
    dirname = 'FAMILY FILES';
end
if exist(dirname,'dir') == 0
    mkdir(dirname)
    cd(dirname)
else
    cd(dirname)
end
%% Checking for family file
filename = strcat('family',int2str(option),'.mat');
if exist(filename,'file')
    columns = (2*ceil(d) + 1)^2;
    load(filename);
    colFile = size(family,2);
    existence = true;
    if colFile ~= columns || length(x) ~= size(family,1)
        disp('Corresponding family file is not corresponding to the current mesh. Erasing it and generating a new one.')
        existence = false;
    else
        cd ../
    end
else
    existence = false;
end
if ~existence
    N = size(x,1);
    lengthx = max(x(:,1)) - min(x(:,1));
    lengthy = max(x(:,2)) - min(x(:,2));
    h = x(2,1) - x(1,1); % Spacing of the grid
    numx = floor(lengthx/h+1/2)+1; % Number of collumns in the mesh
    numy = floor(lengthy/h+1/2)+1; % Number of rows in the mesh
    columns = (2*ceil(d) + 1)^2;
    family = zeros(N,columns); % Instatiate the family matrix. We arbitrarily say that the number of points inside the neighbourhood doesn't surpass 1/5 of the mesh points.
    partialAreas = family; % Also Instatiate familyInfo matrix with zeros
    maxNeigh = 0;
    switch PA_alg
        case "PA-AC"
            %% _____ Algorithm 1: compute interactions_______
            for iI = 1:N
               set = 1;
               % Compute the maximum number of one-sided neighbors interactions on the x-direction
               Nx = floor(horizon/h + 1/2 - 1e-14);
               % Check for borders
               if rem(iI,numx) == 0
                   i = numx;
               else
                   i = rem(iI,numx);
               end  
               if i-1 < Nx
                  Nx_esq = i-1;
                  Nx_dir = Nx;
               elseif numx - i < Nx
                  Nx_esq = Nx;
                  Nx_dir = numx - i;
               else
                  Nx_esq = Nx;
                  Nx_dir = Nx;
               end
               j = ceil(iI/numx); % x = (x_i,y_j); first line, j = 1; last line, j = num
               for k = i-Nx_esq:i+Nx_dir
                   % Compute maximum number of one-sided neighbor interaction along the y-direction Ny
                   if k == i
                      Ny = Nx;
                   else
                      norma_1 = abs(x(iI,1) - x(k + (j-1)*numx,1));
                      Ny = floor((sqrt(horizon^2 - (norma_1 - h/2)^2))/h + 1/2 - 1e-14);
                   end
                   % Check for borders
                   if j-1 < Ny 
                      Ny_bot = j-1;
                      Ny_top = Ny;
                   elseif numy - j < Ny
                      Ny_bot = Ny;
                      Ny_top = numy - j;
                   else
                      Ny_bot = Ny;
                      Ny_top = Ny;
                   end
                   for l = j - Ny_bot: j + Ny_top
                       if (l-1)*numx + k ~= iI
                        family(iI,set) = (l-1)*numx+k;
                        set = set + 1;
                       end
                   end
               end
               if set > maxNeigh
                   maxNeigh = set;
               end
            %% _____ Algorithm 2: compute partial areas___
                   % Map neighbor cell to top-right quadrant
                   % Not necessary; x_tr = zeros(size(x)); % x translated
               N_family = set - 1; % Number of points whithin the neighbourhood of i
               for setII = 1:N_family
                   iII = family(iI,setII);
                   if x(iII,1) < x(iI,1)           
                      x_tr = x(iI,1) + (x(iI,1) - x(iII,1)); % Translate cell from left to right
                   else
                      x_tr = x(iII,1);
                   end
                   if x(iII,2) < x(iI,2)
                      y_tr = x(iI,2) + (x(iI,2) - x(iII,2)); % Translate cell from bottom to top
                   else
                      y_tr = x(iII,2);
                   end
                   if x(iII,1) == x(iI,1)
                      % Rotate cell from +y-axis to +x-axis
                      Delta = y_tr - x(iI,2);
                      y_tr = x(iI,2);
                      x_tr = x(iI,1) + Delta;
                   end
                   % Count number of corners of t_k inside the neighborhood of i
                   counter = 0;
                   for n = -1:2:1
                      for m = -1:2:1
                        if (x_tr + n*h/2 - x(iI,1))^2 + (y_tr + m*h/2 - x(iI,2))^2 < horizon^2
                            counter = counter + 1;
                        end
                      end
                   end
                   % Compute partial area
                   if counter == 4
                      partialAreas(iI,setII) = h^2; % CASE I
                   elseif counter == 3 
                       % Compute A: CASE II
                       % Geometric parameters
                       H1 = y_tr + h/2 - x(iI,2); L1 = sqrt(horizon^2 - H1^2); L2 = x_tr + h/2 - x(iI,1); H2 = sqrt(horizon^2 - L2^2); d = sqrt((H1-H2)^2 + (L2 - L1)^2); l = sqrt(horizon^2 - (d/2)^2); gamma = asin(d/2/horizon);
                       % Area
                       A_k = h^2 - 1/2*(L2-L1)*(H1 - H2) + gamma*horizon^2 - 1/2*d*l;
                       % Family Info
                       partialAreas(iI,setII) = A_k; % CASE II
                   elseif counter == 2
                       if y_tr == x(iI,2)
                               if x_tr + h/2 > x(iI,1) + horizon - 1e-14
                                   % Compute A: CASE III(a2)
                                   % Geometric parameters
                                   l = sqrt(horizon^2 - (h/2)^2); gamma = asin(h/2/horizon);
                                   % Area
                                   A_k = (l - (x_tr - h/2 - x(iI,1)))*h + (gamma*horizon^2 - 1/2*h*l);
                                   % Family Info
                                   partialAreas(iI,setII) = A_k; % CASE III(a2)
                               else
                                   % Compute A: CASE III(b)
                                   % Geometric parameters
                                   l = sqrt(horizon^2 - (h/2)^2); L = x_tr + h/2 - x(iI,1); gamma = acos(l/horizon); beta = acos(L/horizon); d = 2*sqrt(horizon^2 - L^2);
                                   % Area
                                   A_k = h^2-(L-l)*h + ((gamma*horizon^2 - 1/2*h*l) - (beta*horizon^2 - 1/2*d*L));
                                   % Family Info
                                   partialAreas(iI,setII) = A_k; % CASE III(b)
                               end                   
                       else
                               rbr_sq = (x_tr + h/2 - x(iI,1))^2 + (y_tr - h/2 - x(iI,2))^2; % Distance squared of the bottom-right corner of t_k from i
                               if rbr_sq > horizon^2
                                  % Compute A: CASE III(a1)
                                  % Geometric parameters
                                  H1 = y_tr + h/2 - x(iI,2); H2 = y_tr - h/2 - x(iI,2); L1 = sqrt(horizon^2 - H1^2); L2 = sqrt(horizon^2 - H2^2); d = sqrt((L2 - L1)^2 + h^2); l = sqrt(horizon^2 - (d/2)^2); gamma = asin(d/2/horizon);
                                  % Area
                                  A_k = h*((L1+L2)/2 - (x_tr - h/2 - x(iI,1))) + gamma*horizon^2 - 1/2*d*l;
                                  % Family Info
                                   partialAreas(iI,setII) = A_k; % CASE III(a1)
                               else
                                  % Compute A: CASE III(c)
                                  % Geometric parameters
                                  L1 = x_tr - h/2 - x(iI,1); L2 = x_tr + h/2 -x(iI,1); H1 = sqrt(horizon^2 - L1^2); H2 = sqrt(horizon^2 - L2^2); d = sqrt(h^2 + (H1-H2)^2); l = sqrt(horizon^2 - (d/2)^2); gamma = asin(d/2/horizon);
                                  % Area
                                  A_k = h*((H1+H2)/2 - (y_tr - h/2 - x(iI,2))) + gamma*horizon^2 - 1/2*d*l;
                                  % Family Info
                                  partialAreas(iI,setII) = A_k; % CASE III(c)
                               end                   
                       end           
                   elseif counter == 1
                       % Compute A: CASE IV
                       % Geometric parameters
                       L1 = x_tr - h/2 - x(iI,1); H2 = y_tr - h/2 - x(iI,2); H1 = sqrt(horizon^2 - L1^2); L2 = sqrt(horizon^2 - H2^2); d = sqrt((L2-L1)^2 + (H1 - H2)^2); l = sqrt(horizon^2 - (d/2)^2); gamma = asin(d/2/horizon);
                       % Area
                       A_k = 1/2*(L2-L1)*(H1-H2) + gamma*horizon^2 - 1/2*d*l;
                       % Generate family info
                       partialAreas(iI,setII) = A_k; % CASE IV
                   else
                       if x_tr - h/2 < x(iI,1) + horizon
                           %Compute A: CASE V
                           % Geometric parameters
                           l = x_tr - h/2 - x(iI,1); d= 2*sqrt(horizon^2 - l^2); gamma = acos(l/horizon);
                           % Area
                           A_k = gamma*horizon^2 - 1/2*d*l;
                           % Family Info
                           partialAreas(iI,setII) = A_k; % CASE V
                       else
                           partialAreas(iI,setII) = 0;
                       end           
                   end
               end
               display(strcat('Generating family..',num2str(iI/N*100),'%'));
            end
        case "PA-HHB"
            %% Algorithm: HHB
            for iI = 1:N
                x_i = x(iI,:);
                neigh_index = 1;
                for iII = 1:N
                    if iII ~= iI
                        x_j = x(iII,:);
                        xi  = x_j - x_i;
                        norma = norm(xi);
                        % Check if cell tau_k is contained within the neighborhood of i
                        if norma + h/2 <= horizon
                            family(iI,neigh_index) = iII;
                            partialAreas(iI,neigh_index) = h^2;
                            neigh_index = neigh_index+1;
                        elseif norma - h/2 <= horizon
                            family(iI,neigh_index) = iII;
                            partialAreas(iI,neigh_index) = 1/h*(horizon - (norma-h/2))*h^2;
                            neigh_index = neigh_index+1;
                        end
                    end
                end
                display(strcat('Generating family..',num2str(iI/N*100),'%'));
            end
        otherwise
             for iI = 1:N
                x_i = x(iI,:);
                neigh_index = 1;
                for iII = 1:N
                    if iII ~= iI
                        x_j = x(iII,:);
                        xi  = x_j - x_i;
                        norma = norm(xi);
                    if norma <= horizon
                            family(iI,neigh_index) = iII;
                            partialAreas(iI,neigh_index) = h^2;
                            neigh_index = neigh_index+1;
                    end
                    end
                end
             end
%            error("Partial area algorithm not yet implemented")
    end
    %% SAVE
    if find(family(:,columns)~=0)~=0
        disp('Warning: family matrix is not big enough to acomodate all neighbours. Please check the file generateFamily.m.');
    else
        save(filename,'family','partialAreas','maxNeigh')
    end
    cd ../
end
%% Surface correction algorithm
    % Volume method
    surfaceCorrection = zeros(size(partialAreas));
    V0 = pi*horizon^2;
     N = size(x,1);
    for iI = 1:N
        for iII = 1:length(family(iI,family(iI,:)~=0))
            jj = family(iI,iII);
            switch surf_Alg
                case "Volume"
                    surfaceCorrection(iI,iII) = 2*V0/(sum(partialAreas(iI,:)) + sum(partialAreas(jj,:)));
                otherwise
                    surfaceCorrection(iI,iII) = 1;
            end
        end
        %display(strcat('Generating surface correction..',num2str(iI/N*100),'%'));
    end
end

