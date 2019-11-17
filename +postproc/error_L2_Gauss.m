% Function to evaluate the error in the L2 norm
% INPUTS:
% - uh: discretized displacement
% - u : analytical displacement
% - x : x - position of the nodes
% - y : y - position of the nodes

function [error] = error_L2_Gauss(x,y,uh,vh,X,h)
    Nu = size(X,1);
    RIntU = 0;
    RIntV = 0;
    for i = 1:Nu
        if (X(i,1) ~= 0 && X(i,1) ~= 1) && (X(i,2) ~= 0 && X(i,2) ~= 1)
            [IntU,IntV] = gaussQuadrature2D(uh(i),vh(i),h,X(i,1),X(i,2));
            RIntU = RIntU + IntU;
            RIntV = RIntV + IntV;
        end
    end
    error = sqrt(RIntU + RIntV);
end

function [U_Int,V_Int] = gaussQuadrature2D(uh_i,vh_i,h,c1,c2)
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
    [uex(g,l),vex(g,l)] = u_exata(c1+h/2*alpha(g),c2 + h/2*alpha(l));
    end
end
% Computing the integral
U_Int = 0;
V_Int = 0;
for k=1:3
    for j=1:3
        U_Int = U_Int+((uh_i - uex(k,j))^2*A(k)*A(l))*h^2/4;
        V_Int = V_Int+((vh_i - vex(k,j))^2*A(k)*A(l))*h^2/4;
    end
end

end




function [U_Int,V_Int] = gaussQuadratureX(c,uh_i,vh_i,h,y)
% INPUT: 
%c: Center (x_i)
%uh_i: Numerical u displacement on x_i
%vh_i: Numerical v displacement on x_i
%h: interval of the mesh
%y: y position of the node

% Gauss parameters
    alpha1 = -0.7745966692; A1 = 0.555555555;
    alpha2 = 0;             A2 = 0.888888888;
    alpha3 = 0.7745966692;  A3 = 0.555555555;
% Computing analytical values
    [uex1,vex1] = u_exata(c+h/2*alpha1,y);
    [uex2,vex2] = u_exata(c+h/2*alpha2,y);
    [uex3,vex3] = u_exata(c+h/2*alpha3,y);
% Computing the integral
    U_Int = ((uh_i - uex1)^2 * A1 + (uh_i - uex2)^2 * A2 + (uh_i - uex3)^2 * A3)*h/2;
    V_Int = ((vh_i - vex1)^2 * A1 + (vh_i - vex2)^2 * A2 + (vh_i - vex3)^2 * A3)*h/2;
end

function [uex,vex] = u_exata(x,y)

global U20 U02 U11 U10 U01 U00 V20 V02 V11 V10 V01 V00
% U
uex = U20*x^2 + U02 *y^2 + U11*x*y + U10*x + U01 * y + U00;
% V
vex = V20*x^2 + V02 *y^2 + V11*x*y + V10*x + V01 * y + V00;
end

function main_old() % Just to keep old code
% First integral: x - Gauss Quadrature - Gauss Legendre
    RInt1 = [];
    RInt2 = [];
    ynew = [];
    for j = 1:size(y,2)
        if y(j) ~= -1 && y(j) ~= 1 % Bounds
            uh_x = uh(j,:); % Choosing a line in uh (that is a matrix) and assigning it to a vector uh_x
            vh_x = vh(j,:);
            Rint = [0 0];
            for i = 1:size(uh_x,2) % size(uh_x,2) = size(vh_x,2)
                if x(i) ~= -1 && x(i) ~= 1 
                    [IntegralU,IntegralV] = gaussQuadratureX(x(i),uh_x(i),vh_x(i),hx,y(j));
                    Rint = Rint  + [IntegralU,IntegralV];
                end
            end
            RInt1 = [RInt1;Rint];
            ynew =[ynew y(j)];
        end
    end
    RInt1  = RInt1';
    % Second integral: y - Gauss Quadrature - Gauss Legendre
    RInt2_U = trapz(ynew,RInt1(1,:));
    RInt2_V = trapz(ynew,RInt1(2,:));
    errorU = sqrt(RInt2_U);
    errorV = sqrt(RInt2_V);
end
