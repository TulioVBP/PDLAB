% Function to evaluate the error in the L2 norm
% INPUTS:
% - uh: discretized displacement
% - u : analytical displacement
% - x : x - position of the nodes
% - y : y - position of the nodes

function [errorU,errorV] = error_H1_Gauss(x,y,uh,vh,X,h)
    Nu = size(X,1);
    RIntU = 0;
    RIntV = 0;
    M = size(x,2);
    for i = 1:Nu
        if (X(i,1) ~= 0 && X(i,1) ~= 1) && (X(i,2) ~= 0 && X(i,2) ~= 1)  % Only the inner nodes
            [DxU,DyU,DxV,DyV] = grad_u_num(uh,vh,i,M,h)
            [IntU,IntV] = gaussQuadrature2D(uh(i),vh(i),h,X(i,1),X(i,2),DxU,DyU,DxV,DyV);
            RIntU = RIntU + IntU;
            RIntV = RIntV + IntV;
        end
    end
    errorU = sqrt(RIntU);
    errorV = sqrt(RIntV);
end

function [U_Int,V_Int] = gaussQuadrature2D(uh_i,vh_i,h,c1,c2,uh_x,uh_y,vh_x,vh_y)
% INPUT: 
% c1: Center (x_i)
% c2: Center (y_i)
% uh_i: Numerical u displacement on x_i
% vh_i: Numerical v displacement on x_i
% h: interval of the mesh
% y: y position of the node

% Gauss parameters
    alpha(1) = -0.7745966692; A(1) = 0.555555555;
    alpha(2) = 0;             A(2) = 0.888888888;
    alpha(3) = 0.7745966692;  A(3) = 0.555555555;
% Computing analytical values
for g = 1:3
    for l = 1:3 
    [uex(g,l),vex(g,l)] = u_exata(c1+h/2*alpha(g),c2 + h/2*alpha(l)); % f
    [uex_x(g,l),uex_y(g,l),vex_x(g,l),vex_y(g,l)] = grad_u_exata(c1+h/2*alpha(g),c2 + h/2*alpha(l)); % grad F
    end
end
% Computing the integral
U_Int = 0;
V_Int = 0;
for k=1:3
    for j=1:3
        U_Int = U_Int+((uh_i - uex(k,j))^2+(uh_x-uex_x(k,j))^2 + (uh_y-uex_y(k,j))^2)*A(k)*A(l)*h^2/4;
        V_Int = V_Int+((vh_i - vex(k,j))^2 + (vh_x-vex_x(k,j))^2 + (vh_y-vex_y(k,j))^2)*A(k)*A(l)*h^2/4;
    end
end

end

%% Exact solution
function [uex,vex] = u_exata(x,y)

global U20 U02 U11 U10 U01 U00 V20 V02 V11 V10 V01 V00
% U
uex = U20*x^2 + U02 *y^2 + U11*x*y + U10*x + U01 * y + U00;
% V
vex = V20*x^2 + V02 *y^2 + V11*x*y + V10*x + V01 * y + V00;
end

%% Exact gradient
function [DxU,DyU,DxV,DyV] = grad_u_exata(x,y)
global U20 U02 U11 U10 U01 V20 V02 V11 V10 V01
% Ux
DxU = 2*U20*x + U11*y + U10;
% Uy
DyU = 2*U02*y + U11*x + U01;
% Vx
DxV = 2*V20*x + V11*y + V10;
% Vy
DyV = 2*V02*y + V11*x + V01;
end

%% Numerical gradient

function [DxU,DyU,DxV,DyV] = grad_u_num(uh,vh,i,M,h)
    DxU = (uh(i+1) - uh(i-1))/(2*h);
    DyU = (uh(i+M) - uh(i-M))/(2*h);
    DxV = (vh(i+1) - vh(i-1))/(2*h);
    DyV = (vh(i+M) - vh(i-M))/(2*h);
end
