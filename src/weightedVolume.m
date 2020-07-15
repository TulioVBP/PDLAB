 function m_anl = weightedVolume(parameters)
%% Analytical 'm': m is equal inside all the domain, as if it was on the bulk
% NBL MODEL
% INPUT : 
% x - Position of the nodes
% familySet - familySet(i,j) indicates the index of each node in the
%             neighborhood of node 'i'
% familyInfo - familyInfo(i,j) indicates whether the node is on the
%              boundary or not
% V  - V(i) represents the element volume of each node (i)
% parameters = [horizon option alfa]
% ------
% OUTPUT : 
% m_anl - Weighted Volume
% ------
%global alfa
horizon = parameters(1); 
option = parameters(2); 
alfa = parameters(3);
    %N = size(x,1); % Number of nodes
    %m = ones(N,1); % Initializing the "m" vector
    switch option
        case 1
            l = horizon/3;
            if alfa == 0
                m_anl = l^2 * pi *(l^2 - exp(-(horizon^2/l^2))*(l^2 + horizon^2));
            elseif alfa == 1
                m_anl = pi^(3/2)/2*(l)^3;
            end
        case 2
            %P0
            if alfa == 0
                 m_anl = pi*horizon^4/2;
            elseif alfa == 1
                 m_anl = 2*pi*horizon^3/3;
            end
        case 3
            % P1
            if alfa == 0
                 m_anl = pi*horizon^4/10;
            elseif alfa == 1
                 m_anl = pi*horizon^3/6;
            end
        case 4
            % P3
            if alfa == 0
                 m_anl = pi*horizon^4/14;
            elseif alfa == 1
                 m_anl = 2*pi*horizon^3/15;
            end
        case 5
            if alfa == 0
                 m_anl = 5*pi*horizon^4/84;
            elseif alfa == 1
                 m_anl = 5*pi*horizon^3/42;
            end
        case 6
            if alfa == 0
                 m_anl = 7*pi*horizon^4/132;
            elseif alfa == 1
                 m_anl = pi*horizon^3/9;
            end
        case 7
            % Singular
            if alfa == 0
                 m_anl = pi*horizon^4/6;
            elseif alfa == 1
                 m_anl = pi*horizon^3/3;
            end
    end
    %m = m_anl*m;
end
