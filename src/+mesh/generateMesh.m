function [x,A] = generateMesh(h,geometry,type)
% [x,A] = generateMesh(h,geometry,type) is a function that takes a the grid
% spacing, the geometry vector [height, variables] and the type of the mesh
% ('inclusive', expansive')
if nargin < 3
    type = 'expansive';
end
a = geometry(1); % height [m]
b = geometry(2); % length [m]
switch type
    case 'expansive'
        % Retangular
        N = floor(a/h+1/2)+1; % Number of rows
        M = floor(b/h+1/2)+1; % Number of collumns
        x = zeros(N*M,2);
        for ii = 1:N
            for jj = 1:M
                x((ii-1)*M+jj,:) = [(jj-1)*h,(ii-1)*h];
            end
        end
    case 'inclusive'
        % Retangular
        N = floor(a/h+1/2); % Number of rows
        M = floor(b/h+1/2); % Number of collumns
        x = zeros(N*M,2);
        for ii = 1:N
            for jj = 1:M
                x((ii-1)*M+jj,:) = [(jj-1)*h+h/2,(ii-1)*h+h/2];
            end
        end
    case default 
    disp("Type "+type+ " doesn't exist.")
end

A = h^2; % Elements' area
% Symmetry
x = x - [(max(x(:,1))+min(x(:,1)))/2 (max(x(:,2))+min(x(:,2)))/2];

end