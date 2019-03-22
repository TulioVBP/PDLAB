function m = weightedVolume(x,horizon,option)
%% Analytical 'm': m is equal inside all the domain, as if it was on the bulk
% NBL MODEL
% INPUT : 
% x - Position of the nodes
% familySet - familySet(i,j) indicates the index of each node in the
%             neighborhood of node 'i'
% familyInfo - familyInfo(i,j) indicates whether the node is on the
%              boundary or not
% V  - V(i) represents the element volume of each node (i)
% horizon - Peridynamic horizon
% option - defines which influence function is used
% ------
% OUTPUT : 
% m - Weighted Volume
% ------
global alpha
    N = size(x,1); % Number of nodes
    m = ones(N,1); % Initializing the "m" vector
    switch option
        case 1
            l = horizon/4.3;
            m_anl = l^2 * pi *(l^2 - exp(-(horizon^2/l^2))*(l^2 + horizon^2));
        case 2
            %P0
            if alpha == 0
                 m_anl = pi*horizon^4/2;
            elseif alpha == 1
                 m_anl = 2*pi*horizon^3/3;
            end
        case 3
            % P1
            if alpha == 0
                 m_anl = pi*horizon^4/10;
            elseif alpha == 1
                 m_anl = pi*horizon^3/6;
            end
        case 4
            % P3
            if alpha == 0
                 m_anl = pi*horizon^4/14;
            elseif alpha == 1
                 m_anl = 2*pi*horizon^3/15;
            end
        case 5
            if alpha == 0
                 m_anl = 5*pi*horizon^4/84;
            elseif alpha == 1
                 m_anl = 5*pi*horizon^3/42;
            end
        case 6
            if alpha == 0
                 m_anl = 7*pi*horizon^4/132;
            elseif alpha == 1
                 m_anl = pi*horizon^3/9;
            end
    end
    m = m_anl*m;
end
