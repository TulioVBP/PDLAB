function [x,A] = generateMesh(h,geometry)
a = geometry(1); % height [m]
b = geometry(2); % length [m]
% Retangular
N = floor(a/h+1/2) + 1; % Number of rows
M = floor(b/h+1/2) + 1; % Number of collumns
x = zeros(N*M,2);
for ii = 1:N
    for jj = 1:M
        x((ii-1)*M+jj,:) = [(jj-1)*h,(ii-1)*h];
    end
end
A = h^2; % Elements' area
end