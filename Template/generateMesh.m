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
% Hole
xc = [20,20;
      36.5,51;
      20,100]*10^-3;
rc = [5,10,5]*10^-3;
for ll = 1:size(xc,1)
    r = vecnorm((x-xc(ll,:))')';
    hole = r < rc(ll);
    x = x(~hole,:);
end

A = h^2; % Elements' area
% Symmetry
x = x - [max(x(:,1))/2 max(x(:,2))/2];
end