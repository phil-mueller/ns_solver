
function [stecil,b] = stamp(i, j, space, FV, boundary, value, str) 
%stecil calculate the linear equation for node (i, j)
%
%  input:
%    i         node number in space.X direction
%    j         node number in space.Y direction
%    space.X         space.X position of the nodes
%    space.Y         space.Y position of the nodes
%    b         right-hand side value for node (i,j)
%    boundary  defines the boundary conditions
%
%  output:
%    stecil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init
b=0;
n = size(space.X, 1);
m = size(space.X, 2);
stecil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon
nodePosition = 'inner Node';

if i==n
    nodePosition = 'South';
%     ib.u(i,j) = false;
%     ib.v(i,j) = false;
elseif i==1
    nodePosition = 'North';
%     ib.u(i,j) = false;
%     ib.v(i,j) = false;
end

if j == 1
    nodePosition = 'West';
%     ib.u(i,j) = false;
%     ib.v(i,j) = false;
    if i == 1
        nodePosition = 'WestNorth';
%         ib.u(i,j) = false;
%         ib.v(i,j) = false;
    elseif i == n
        nodePosition = 'WestSouth';
%         ib.u(i,j) = false;
%         ib.v(i,j) = false;
    end
    
elseif j == m
    nodePosition = 'East';
%     ib.p(i,j) = false;
    
    if i == 1
        nodePosition = 'EastNorth';
%         ib.p(i,j) = false;
    elseif i == n
        nodePosition = 'EastSouth';
%         ib.p(i,j) = false;
    end
end

% If an obstacle with corners shall be used, define new boundary condition
% points

% Calculate the equation for the correct node position
switch nodePosition
    
    case 'inner Node'
        
        
        % Nomenclature:
        %
        %    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)
        %
        %                 |                 |
        %
        %       nW - - - - nw ------ n ------ ne - - - nE
        %                 |                 |
        %       |         |        |        |       |
        %                 |                 |
        %   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
        %                 |                 |
        %       |         |        |        |       |
        %                 |                 |
        %      sW - - - - sw ------ s ------ se - - - sE
        %
        %                 |                 |
        %
        %   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
        %
        % Indexing of stecil:
        
        %    D_4 - D_1 - D2
        %     |     |     |
        %    D_3 - D_0 - D3
        %     |     |     |
        %    D_2 -  D1 - D4
        
        
        
        y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
        y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
        y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
        y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
        y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
        y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
        y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
        y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
        y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
        
        y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
        y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
        y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
        y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
        
        y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
        y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
        y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
        y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
        y_nE =  (y_NE + y_E)/2;  x_nE =  (x_NE + x_E)/2;
        y_sE =  (y_E +  y_SE)/2; x_sE =  (x_E +  x_SE)/2;
        
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
        y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
        y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
        
        
        
        % Around s
        
        dy_Sw_Se = y_Se - y_Sw;  dx_Sw_Se = x_Se - x_Sw;
        dy_Se_e = y_e - y_Se;    dx_Se_e = x_e - x_Se;
        dy_e_w = y_w - y_e;      dx_e_w = x_w - x_e;
        dy_w_Sw = y_Sw - y_w;    dx_w_Sw = x_Sw - x_w;
        
        % Around e
        
        dy_s_sE = y_sE - y_s;    dx_s_sE = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;  dx_sE_nE = x_nE - x_sE;
        dy_nE_n = y_n - y_nE;    dx_nE_n = x_n - x_nE;
        dy_n_s = y_s - y_n;      dx_n_s = x_s - x_n;
        
        % Around n
        
        dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;
        dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;  dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w - y_Nw;   dx_Nw_w  = x_w - x_Nw;
        
        
        % Around w
        dy_sW_s = y_s - y_sW;    dx_sW_s = x_s - x_sW;
        dy_s_n  = y_n - y_s;     dx_s_n  = x_n - x_s;
        dy_n_nW = y_nW - y_n;    dx_n_nW = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;  dx_nW_sW = x_sW - x_nW;
        
        
        % Around P
        
        dy_sw_se = y_se - y_sw;   dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;   dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;   dx_nw_sw = x_sw - x_nw;
        
        % Areas
        
        
        
        S_P                 = abs( x_ne*y_se - y_ne*x_se...
            + x_se*y_sw - y_se*x_sw...
            + x_sw*y_nw - y_sw*x_nw...
            + x_nw*y_ne - y_nw*x_ne) / 2;
        S_s                 = abs( x_e*y_Se - y_e*x_Se...
            + x_Se*y_Sw - y_Se*x_Sw...
            + x_Sw*y_w - y_Sw*x_w...
            + x_w*y_e - y_w*x_e) / 2;
        S_e                 = abs( x_nE*y_sE - y_nE*x_sE...
            + x_sE*y_s - y_sE*x_s...
            + x_s*y_n - y_s*x_n...
            + x_n*y_nE - y_n*x_nE) / 2;
        S_n                 = abs( x_Ne*y_e - y_Ne*x_e...
            + x_e*y_w - y_e*x_w...
            + x_w*y_Nw - y_w*x_Nw...
            + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
        S_w                 = abs( x_n*y_s - y_n*x_s...
            + x_s*y_sW - y_s*x_sW...
            + x_sW*y_nW - y_sW*x_nW...
            + x_nW*y_n - y_nW*x_n) / 2;
        
        
        
        %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
        
        build_inner
        
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        % East
        stecil(index(i, j+1)) = D3;
        % West
        stecil(index(i, j-1)) = D_3;
        % South
        stecil(index(i+1, j)) = D1;
        % North
        stecil(index(i-1, j)) = D_1;
        % NW
        stecil(index(i-1, j-1)) = D_4;
        % NE
        stecil(index(i-1, j+1)) = D2;
        % SW
        stecil(index(i+1, j-1)) = D_2;
        % SE
        stecil(index(i+1, j+1)) = D4;
        % P
        stecil(index(i, j)) = D0;
        
    case 'South'
%         if strcmp(boundary.south, 'Neumann')
            
            % Nomencature:
            %    NW=(i-1,j-1)  Nw - N=(i-1,j)- Ne    NE=(i-1,j+1)
            %
            %                  |                |
            %
            %       nW- - - -  nw----- n ------ne - - - nE
            %                  |                |
            %      etaW        etaw   eta    etae      etaE
            %                  |                |
            %    W=(i,j-1) - - w - - P=(i,j) - -e - - E=(i,j+1)
            
            
            
            
            y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
            y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
            y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
            y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
            y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
            y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
            
            
            y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
            y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
            
            
            y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_nE =  (y_NE + y_E)/2;  x_nE =  (x_NE + x_E)/2;
            
            
            y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
            y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
            
            y_etaW = (y_nW + y_W)/2; x_etaW = (x_nW + x_W)/2;
            y_etaw = (y_nw + y_w)/2; x_etaw = (x_nw + x_w)/2;
            y_eta  = (y_n  + y_P)/2; x_eta  = (x_n  + x_P)/2;
            y_etae = (y_ne + y_e)/2; x_etae = (x_ne + x_e)/2;
            y_etaE = (y_nE + y_E)/2; x_etaE = (x_nE + x_E)/2;
            
            
            % Around etaw
            
            dy_W_P  = y_P - y_W;    dx_W_P  = x_P - x_W;
            dy_P_n  = y_n - y_P;    dx_P_n  = x_n - x_P;
            dy_n_nW = y_nW - y_n;   dx_n_nW = x_nW - x_n;
            dy_nW_W = y_W - y_nW;   dx_nW_W = x_W - x_nW;
            
            
            % Around etae
            
            dy_P_E  = y_E - y_P;    dx_P_E  = x_E - x_P;
            dy_E_nE = y_nE - y_E;   dx_E_nE = x_nE - x_E;
            dy_nE_n = y_n - y_nE;   dx_nE_n = x_n - x_nE;
            dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
            
            % Around n
            
            dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;
            dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
            dy_Ne_Nw = y_Nw - y_Ne;  dx_Ne_Nw = x_Nw - x_Ne;
            dy_Nw_w  = y_w - y_Nw;   dx_Nw_w  = x_w - x_Nw;
            
            
            % Around P
            
            dy_e_ne = y_ne - y_e;     dx_e_ne = x_ne - x_e;
            dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
            dy_nw_w =  y_w - y_nw;    dx_nw_w =  x_w - x_nw;
            
            
            % Areas
            
            
            S_eta                 = abs( x_ne*y_e - y_ne*x_e...
                + x_e*y_w - y_e*x_w...
                + x_w*y_nw - y_w*x_nw...
                + x_nw*y_ne - y_nw*x_ne) / 2;
            S_etae                 = abs( x_nE*y_E - y_nE*x_E...
                + x_E*y_P - y_E*x_P...
                + x_P*y_n - y_P*x_n...
                + x_n*y_nE - y_n*x_nE) / 2;
            S_n                 = abs( x_Ne*y_e - y_Ne*x_e...
                + x_e*y_w - y_e*x_w...
                + x_w*y_Nw - y_w*x_Nw...
                + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
            S_etaw                 = abs( x_n*y_P - y_n*x_P...
                + x_P*y_W - y_P*x_W...
                + x_W*y_nW - y_W*x_nW...
                + x_nW*y_n - y_nW*x_n) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_south
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % West
            stecil(index(i, j-1)) = D_3;
            % North
            stecil(index(i-1, j)) = D_1;
            % NW
            stecil(index(i-1, j-1)) = D_4;
            % NE
            stecil(index(i-1, j+1)) = D2;
            % P
            stecil(index(i, j)) = D0;
            
            if strcmp(str,'p')
                if strcmp(boundary.p.south, 'neumann') || strcmp(boundary.p.south, 'wall')
                    b = 0;
                elseif strcmp(boundary.p.south, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.p_out;
                end
            elseif strcmp(str, 'u')
                if strcmp(boundary.u.south, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.u.south, 'neumann')
                    b = 0;
                elseif strcmp(boundary.u.south, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.u_in;
                end
            elseif strcmp(str, 'v')
                if strcmp(boundary.v.south, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.v.south, 'neumann')
                    b = 0;
                elseif  strcmp(boundary.v.south, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.v_in;
                end
            else
                disp('Enter a valid string');
            end
%         elseif strcmp(boundary.u.south, 'wall') || strcmp(boundary.v.south, 'wall')
%             stecil(index(i, j)) = 1;
%             b = 0.0;
%         else
%             warning('boundary condition not implemented')
%         end
        
    case 'North'
        
%         if strcmp(boundary.u.north, 'Robin') || strcmp(boundary.v.north, 'Robin') || strcmp(boundary.u.north, 'Neumann') ||...
%                 strcmp(boundary.v.north, 'Neumann')
            

            
            % Nomencature:
            %    W=(i,j-1) --  w - -P=(i,j) - - e - - E=(i,j+1)
            %                  |                |
            %   sigmaW    sigmaw    sigma    sigmae    sigmaE
            %                  |                |
            %      sW - - - -  sw----- s ------se - - - sE
            %
            %                  |                |
            %
            %   SW=(i+1,j-1)   Sw  S=(i+1,j) - Se     SE=(i+1,j+1)
            
            
            % Calculate some help values
            
            
            y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
            y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
            y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
            y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
            y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
            y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
            
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
            y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
            y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
            
            y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            y_sE =  (y_E +  y_SE)/2; x_sE =  (x_E +  x_SE)/2;
            
            y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
            y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
            
            
            y_sigmaW = (y_W + y_sW)/2; x_sigmaW = (y_W + y_sW)/2;
            y_sigmaw = (y_w + y_sw)/2; x_sigmaw = (y_w + y_sw)/2;
            y_sigma  = (y_P + y_s)/2;  x_sigma  = (y_P + y_s)/2;
            y_sigmae = (y_e + y_se)/2; x_sigmae = (y_e + y_se)/2;
            y_sigmaE = (y_E + y_sE)/2; x_sigmaE = (y_E + y_sE)/2;
            
            
            
            % Around s
            
            dy_Sw_Se = y_Se - y_Sw; dx_Sw_Se = x_Se - x_Sw;
            dy_Se_e  = y_e - y_Se ; dx_Se_e  = x_e - x_Se ;
            dy_e_w   = y_w - y_e  ; dx_e_w   = x_w - x_e  ;
            dy_w_Sw  = y_Sw - y_w ; dx_w_Sw  = x_Sw - x_w ;
            
            % Around sigmae
            
            dy_s_sE  = y_sE - y_s ; dx_s_sE  = x_sE - x_s ;
            dy_sE_E  = y_E - y_sE ; dx_sE_E  = x_E - x_sE ;
            dy_E_P   = y_P - y_E ;  dx_E_P   = x_P - x_E ;
            dy_P_s   = y_s - y_P ;  dx_P_s   = x_s - x_P ;
            
            % Around sigmaw
            
            dy_sW_s  = y_s - y_sW;  dx_sW_s  = x_s - x_sW;
            dy_s_P   = y_P - y_s ;  dx_s_P   = x_P - x_s ;
            dy_P_W   = y_W - y_P ;  dx_P_W   = x_W - x_P ;
            dy_W_sW  = y_sW - y_W;  dx_W_sW  = x_sW - x_W;
            
            % Around P
            
            dy_sw_se  = y_se - y_sw;  dx_sw_se  = x_se - x_sw;
            dy_se_e   = y_e - y_se ;  dx_se_e   = x_e - x_se ;
            dy_e_w    = y_w - y_e;    dx_e_w    = x_w - x_e;
            dy_w_sw   = y_sw - y_w;   dx_w_sw   = x_sw - x_w;
            
            
            % At boundary
            dl_w_e = norm([dx_e_w; dy_e_w]);
            
            
            
            % Areas
            
            S_sigma                 = abs( x_e*y_se - y_e*x_se...
                + x_se*y_sw - y_se*x_sw...
                + x_sw*y_w - y_sw*x_w...
                + x_w*y_e - y_w*x_e) / 2;
            S_s                 = abs( x_e*y_Se - y_e*x_Se...
                + x_Se*y_Sw - y_Se*x_Sw...
                + x_Sw*y_w - y_Sw*x_w...
                + x_w*y_e - y_w*x_e) / 2;
            S_sigmae                 = abs( x_E*y_sE - y_E*x_sE...
                + x_sE*y_s - y_sE*x_s...
                + x_s*y_P - y_s*x_P...
                + x_P*y_E - y_P*x_E) / 2;
            S_sigmaw                 = abs( x_P*y_s - y_P*x_s...
                + x_s*y_sW - y_s*x_sW...
                + x_sW*y_W - y_sW*x_W...
                + x_W*y_P - y_W*x_P) / 2;
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_north
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % West
            stecil(index(i, j-1)) = D_3;
            % South
            stecil(index(i+1, j)) = D1;
            % SW
            stecil(index(i+1, j-1)) = D_2;
            % SE
            stecil(index(i+1, j+1)) = D4;
            % P
            stecil(index(i, j)) = D0;
            if strcmp(str,'p')
                if strcmp(boundary.p.north, 'neumann') || strcmp(boundary.p.north, 'wall')
                    b = 0;
                elseif strcmp(boundary.p.north, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.p_out;
                end
            elseif strcmp(str, 'u') 
                if strcmp(boundary.u.north, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.u.north, 'neumann')
                    b = 0;
                elseif strcmp(boundary.u.north, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.u_in;
                end
            elseif strcmp(str, 'v')
                if strcmp(boundary.v.north, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.v.north, 'neumann')
                    b = 0;
                elseif strcmp(boundary.v.north, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.v_in;
                end
            else
                disp('Enter a valid string');
            end
                
%             if strcmp(boundary.u.north, 'Neumann') || strcmp(boundary.v.north, 'Neumann')
                % P
%                 stecil(index(i, j)) = D0;
%                 b = q.north*dl_w_e/S_sigma;
%             else
%                 stecil(index(i, j)) = D0 - alpha*dl_w_e/S_sigma;
%                 b = -alpha*Tinf*dl_w_e/S_sigma;
%             end
%             %fbm................................................
%         elseif strcmp(boundary.u.north, 'wall') || strcmp(boundary.v.north, 'wall')
%             stecil(index(i,j)) = 1;
%             b=0.0;
%         else
%             warning('boundary condition not implemented')
%         end
        
        
    case 'East'
        
%         if strcmp(boundary.u.east, 'Robin') || strcmp(boundary.v.east, 'Robin') || strcmp(boundary.u.east, 'Neumann') ||...
%                 strcmp(boundary.v.east, 'Neumann')
            

            
            % Nomencature:
            %   NW=(i-1,j-1)    Nw -Nomega - N=(i-1,j)
            %
            %                   |            |
            %
            %       nW   - - - -nw- nomega - n
            %                   |
            %       |           |            |
            %                   |
            %    W=(i,j-1) - -  w - omega   P=(i,j)
            %                   |
            %       |           |            |
            %                   |
            %       sW - - -    sw- somega - s
            %
            %                   |            |
            %
            %    SW=(i+1,j-1)   Sw -Somega S=(i+1,j)
            
            % Calculate some help values
            
            
            y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
            y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
            y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
            y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
            y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
            y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
            
            y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
            y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
            y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
            
            y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
            y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            
            y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
            y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
            
            
            y_Nomega  = (y_N + y_Nw)/2;  x_Nomega  = (x_N + x_Nw)/2;
            y_nomega  = (y_n  + y_nw)/2; x_nomega  = (x_n  + x_nw)/2;
            y_omega   = (y_P + y_w)/2;   x_omega   = (x_P + x_w)/2;
            y_somega  = (y_s + y_sw)/2;  x_somega  = (x_s + x_sw)/2;
            y_Somega  = (y_S + y_Sw)/2;  x_Somega  = (x_S + x_Sw)/2;
            
            
            % Around somega
            
            dy_Sw_S = y_S - y_Sw;  dx_Sw_S = x_S - x_Sw;
            dy_S_P  = y_P - y_S;   dx_S_P  = x_P - x_S;
            dy_P_w  = y_w - y_P;   dx_P_w  = x_w - x_P;
            dy_w_Sw = y_Sw - y_w;  dx_w_Sw = x_Sw - x_w;
            
            % Around nomega
            
            dy_w_P  = y_P - y_w;   dx_w_P  = x_P - x_w;
            dy_P_N  = y_N - y_P;   dx_P_N  = x_N - x_P;
            dy_N_Nw = y_Nw - y_N;  dx_N_Nw = x_Nw - x_N;
            dy_Nw_w = y_w - y_Nw;  dx_Nw_w = x_w - x_Nw;
            
            % Around w
            
            dy_sW_s  = y_s - y_sW;   dx_sW_s  = x_s - x_sW;
            dy_s_n   = y_n - y_s;    dx_s_n   = x_n - x_s;
            dy_n_nW  = y_nW - y_n;   dx_n_nW  = x_nW - x_n;
            dy_nW_sW = y_sW - y_nW;  dx_nW_sW = x_sW - x_nW;
            
            % Around P
            
            dy_sw_s  = y_s - y_sw;   dx_sw_s  = x_s - x_sw;
            dy_s_n   = y_n - y_s;    dx_s_n   = x_n - x_s;
            dy_n_nw  = y_nw - y_n;   dx_n_nw  = x_nw - x_n;
            dy_nw_sw = y_sw - y_nw;  dx_nw_sw = x_sw - x_nw;
            
            % At boundary
            
            dl_s_n = norm([dx_s_n; dy_s_n]);
            
            % Areas
            
            
            S_omega                 = abs( x_n*y_s - y_n*x_s...
                + x_s*y_sw - y_s*x_sw...
                + x_sw*y_nw - y_sw*x_nw...
                + x_nw*y_n - y_nw*x_n) / 2;
            S_somega                 = abs( x_P*y_S - y_P*x_S...
                + x_S*y_Sw - y_S*x_Sw...
                + x_Sw*y_w - y_Sw*x_w...
                + x_w*y_P - y_w*x_P) / 2;
            
            S_nomega                 = abs( x_N*y_P - y_N*x_P...
                + x_P*y_w - y_P*x_w...
                + x_w*y_Nw - y_w*x_Nw...
                + x_Nw*y_N - y_Nw*x_N ) / 2;
            S_w                 = abs( x_n*y_s - y_n*x_s...
                + x_s*y_sW - y_s*x_sW...
                + x_sW*y_nW - y_sW*x_nW...
                + x_nW*y_n - y_nW*x_n) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_east
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            % West
            stecil(index(i, j-1)) = D_3;
            % South
            stecil(index(i+1, j)) = D1;
            % North
            stecil(index(i-1, j)) = D_1;
            % NW
            stecil(index(i-1, j-1)) = D_4;
            % SW
            stecil(index(i+1, j-1)) = D_2;
            % P
            stecil(index(i, j)) = D0;
            %fbm....................................................
            if strcmp(str,'p')
                if strcmp(boundary.p.east, 'neumann') || strcmp(boundary.p.east, 'wall')
                    b = 0;
                elseif strcmp(boundary.p.east, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.p_out;
                end
            elseif strcmp(str, 'u')
                if strcmp(boundary.u.east, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.u.east, 'neumann')
                    b = 0;
                elseif strcmp(boundary.u.east, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.u_in;
                end
            elseif strcmp(str, 'v')
                if strcmp(boundary.v.east, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.v.east, 'neumann')
                    b = 0;
                elseif strcmp(boundary.v.east, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.v_in;
                end
            else
                disp('Enter a valid string');
            end
%             if strcmp(boundary.u.east, 'Neumann') || strcmp(boundary.v.east, 'Neumann')
%                 % P
%                 stecil(index(i, j)) = D0;
%                 b = q.east*dl_s_n/S_omega; %Equivalent value for u and v??
%             else
%                 stecil(index(i, j)) = D0 - alpha*dl_s_n/S_omega;
%                 b = -alpha*Tinf*dl_s_n/S_omega;
%             end
%             %fbm.....................................................
%         elseif strcmp(boundary.u.east, 'Dirichlet')
%             stecil(index(i,j)) = 1;
%             b=TD.east;
%         else
%             warning('boundary condition not implemented')
%         end
        
        
    case 'West'
        
%         if strcmp(boundary.west, 'Robin') || strcmp(boundary.west, 'Neumann')
            

            
            % Calculate some help values
            
            
            y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
            y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
            y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
            y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
            y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
            y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
            
            y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
            y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
            y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
            
            y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
            y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
            y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
            y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
            
            y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
            y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
            
            
            y_Nepsilon  = (y_N + y_Ne)/2;  x_Nepsilon  = (x_N + x_Ne)/2;
            y_nepsilon  = (y_n  + y_ne)/2; x_nepsilon  = (x_n  + x_ne)/2;
            y_epsilon   = (y_P + y_e)/2;   x_epsilon   = (x_P + x_e)/2;
            y_sepsilon  = (y_s + y_se)/2;  x_sepsilon  = (x_s + x_se)/2;
            y_Sepsilon  = (y_S + y_Se)/2;  x_Sepsilon  = (x_S + x_Se)/2;
            
            
            % Around sepsilon
            
            dy_S_Se  = y_Se - y_S;   dx_S_Se = x_Se - x_S;
            dy_Se_e  = y_e - y_Se;   dx_Se_e = x_e - x_Se;
            dy_e_P   = y_P - y_e;    dx_e_P  = x_P - x_e;
            dy_P_S   = y_S - y_P;    dx_P_S  = x_S - x_P;
            
            % Around nepsilon
            
            dy_P_e   = y_e - y_P;    dx_P_e   = x_e - x_P;
            dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
            dy_Ne_N  = y_N - y_Ne;   dx_Ne_N  = x_N - x_Ne;
            dy_N_P   = y_P - y_N;    dx_N_P   = x_P - x_N;
            
            % Around e
            
            dy_s_sE  = y_sE - y_s;   dx_s_sE  = x_sE - x_s;
            dy_sE_nE = y_nE - y_sE;  dx_sE_nE = x_nE - x_sE;
            dy_nE_n  = y_n - y_nE;   dx_nE_n  = x_n - x_nE;
            dy_n_s   = y_s - y_n;    dx_n_s   = x_s - x_n;
            
            % Around P
            
            dy_s_se  = y_se - y_s;   dx_s_se  = x_se - x_s;
            dy_se_ne = y_ne - y_se;  dx_se_ne = x_ne - x_se;
            dy_ne_n  = y_n - y_ne;   dx_ne_n  = x_n - x_ne;
            dy_n_s   = y_s - y_n;    dx_n_s   = x_s - x_n;
            
            % At boundary
            
            dl_s_n = norm([dx_n_s; dy_n_s]);
            
            % Areas
            
            
            S_epsilon                 = abs( x_s*y_n - y_s*x_n...
                + x_n*y_ne - y_n*x_ne...
                + x_ne*y_se - y_ne*x_se...
                + x_se*y_s - y_se*x_s) / 2;
            
            S_sepsilon                 = abs( x_S*y_P - y_S*x_P...
                + x_P*y_e - y_P*x_e...
                + x_e*y_Se - y_e*x_Se...
                + x_Se*y_S - y_Se*x_S) / 2;
            
            S_nepsilon                 = abs( x_P*y_N - y_P*x_N...
                + x_N*y_Ne - y_N*x_Ne...
                + x_Ne*y_e - y_Ne*x_e...
                + x_e*y_P - y_e*x_P ) / 2;
            
            S_e                 = abs( x_s*y_n - y_s*x_n...
                + x_n*y_nE - y_n*x_nE...
                + x_nE*y_sE - y_nE*x_sE...
                + x_sE*y_s - y_sE*x_s) / 2;
            
            
            %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
            
            build_west
            
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % East
            stecil(index(i, j+1)) = D3;
            % South
            stecil(index(i+1, j)) = D1;
            % North
            stecil(index(i-1, j)) = D_1;
            % NE
            stecil(index(i-1, j+1)) = D2;
            % SE
            stecil(index(i+1, j+1)) = D4;
            % P
            stecil(index(i, j)) = D0;
            %fbm.....................................................
            if strcmp(str,'p')
                if strcmp(boundary.p.west, 'neumann') || strcmp(boundary.p.west, 'wall')
                    b = 0;
                elseif strcmp(boundary.p.west, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.p_out;
                end
            elseif strcmp(str, 'u')
                if strcmp(boundary.u.west, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.u.west, 'neumann')
                    b = 0;
                elseif strcmp(boundary.u.west, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.u_in;
                end
            elseif strcmp(str, 'v')
                if strcmp(boundary.v.west, 'wall')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = 0;
                elseif strcmp(boundary.v.west, 'neumann')
                    b = 0;
                elseif strcmp(boundary.v.west, 'dirichlet')
                    stecil = zeros(1, n*m);
                    stecil(index(i, j)) = 1;
                    b = value.v_in;    
                end
            else
                disp('Enter a valid string');
            end    
            
%             if strcmp(boundary.west, 'Neumann')
%                 % P
%                 stecil(index(i, j)) = D0;
%                 b = q.west*dl_s_n/S_epsilon;
%             else
%                 stecil(index(i, j)) = D0 - alpha*dl_s_n/S_epsilon;
%                 b = -alpha*Tinf*dl_s_n/S_epsilon;
%             end
%             %fbm.....................................................
%         elseif strcmp(boundary.west, 'Dirichlet')
%             stecil(index(i,j)) = 1;
%             b=TD.west;
%         else
%             warning('boundary condition not implemented')
%         end
%         
%         
    case 'EastSouth'
        
%          if strcmp(boundary.south, 'Neumann')
            
%              if strcmp(boundary.east, 'Neumann') || strcmp(boundary.east, 'Robin')
                
                % Nomencature:
                %    NW=(i-1,j-1)  Nw - Nomega  - N=(i-1,j)
                %                  |       |        |
                %                  |       |        |
                %                  |       |        |
                %       nW- - - -  nw----nomega  ----- n
                %                  |       |        |
                %      etaW      etaw   etaomega   eta
                %                  |       |        |
                %    W=(i,j-1) - - w - - omega -- P=(i,j)
                
                
                % Calculate some help values
                
                y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
                y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
                y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
                y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
                
                y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
                y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
                
                y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                
                y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
                
                % Around nomega
                
                dy_w_P  = y_P - y_w;   dx_w_P  = x_P - x_w;
                dy_P_N  = y_N - y_P;   dx_P_N  = x_N - x_P;
                dy_N_Nw = y_Nw - y_N;  dx_N_Nw = x_Nw - x_N;
                dy_Nw_w = y_w - y_Nw;  dx_Nw_w = x_w - x_Nw;
                
                % Around etaw
                
                dy_W_P  = y_P - y_W;    dx_W_P  = x_P - x_W;
                dy_P_n  = y_n - y_P;    dx_P_n  = x_n - x_P;
                dy_n_nW = y_nW - y_n;   dx_n_nW = x_nW - x_n;
                dy_nW_W = y_W - y_nW;   dx_nW_W = x_W - x_nW;
                
                % Around P
                
                dy_n_nw = y_nw - y_n;  dx_n_nw = x_nw - x_n;
                dy_nw_w = y_w - y_nw;  dx_nw_w = x_w - x_nw;
                
                % At boundary
                
                dl_P_n = norm([dx_P_n; dy_P_n]);
                
                % Areas
                
                S_etaw              = abs( x_n*y_P - y_n*x_P...
                    + x_P*y_W - y_P*x_W...
                    + x_W*y_nW - y_W*x_nW...
                    + x_nW*y_n - y_nW*x_n) / 2;
                S_nomega            = abs( x_N*y_P - y_N*x_P...
                    + x_P*y_w - y_P*x_w...
                    + x_w*y_Nw - y_w*x_Nw...
                    + x_Nw*y_N - y_Nw*x_N ) / 2;
                S_etaomega         = abs( x_n*y_P - y_n*x_P...
                    + x_P*y_w - y_P*x_w...
                    + x_w*y_nw - y_w*x_nw...
                    + x_nw*y_n - y_nw*x_n) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_southEast
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                % West
                stecil(index(i, j-1)) = D_3;
                % North
                stecil(index(i-1, j)) = D_1;
                % NW
                stecil(index(i-1, j-1)) = D_4;
                % P
                stecil(index(i, j)) = D0;
%             end
                %fbm.....................................................
                if strcmp(str,'p')
                    if strcmp(boundary.p.south, 'neumann') || strcmp(boundary.p.south, 'wall')
                        if strcmp(boundary.p.east, 'neumann') || strcmp(boundary.p.east, 'wall')
                            b = 0;
                        elseif strcmp(boundary.p.east, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.p_out;  
                        end
                    elseif strcmp(boundary.p.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.p_out;
                    end
                elseif strcmp(str, 'u')
                    if strcmp(boundary.u.south, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.u.south, 'neumann')
                        if strcmp(boundary.u.east, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.u.east, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.u.east, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.u_out;
                        end
                    elseif strcmp(boundary.u.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.u_in;
                    end
                elseif strcmp(str, 'v')
                    if strcmp(boundary.v.south, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.v.south, 'neumann')
                        if strcmp(boundary.v.east, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.v.east, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.v.east, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.v_in;
                        end
                    elseif strcmp(boundary.v.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.v_in;
                    end    
                end                                             
%                 end   
%             if strcmp(boundary.east, 'Neumann')
%                 % P
%                 stecil(index(i, j)) = D0;
%                 b=q.east*dl_P_n/S_etaomega;
%                 
%             elseif strcmp(boundary.east, 'Robin')
%                 % P
%                 stecil(index(i, j)) = D0 - alpha*dl_P_n/S_etaomega;
%                 b = -alpha*Tinf*dl_P_n/S_etaomega;
%                 
%             elseif strcmp(boundary.east, 'Dirichlet')
%               
%                 stecil(index(i, j)) = 1;
%                 b=TD.east;
%             end
% %             end
%             %fbm......................................................
%             
%             
%         elseif strcmp(boundary.south, 'Dirichlet')
%             stecil(index(i, j)) = 1;
%             b=TD.south;
%         end     
    case 'EastNorth'
        
        
%         if strcmp(boundary.north, 'Neumann') || strcmp(boundary.north, 'Robin')
%             
%             if strcmp(boundary.east, 'Neumann') || strcmp(boundary.east, 'Robin')
                
                % Nomencature:
                %
                %  W=(i,j-1) -  -  w---- omega  ---P=(i,j)
                %                  |       |        |
                % sigmaW      sigmaw  sigmaomega   sigma
                %                  |       |        |
                %     sW - - - - -sw - -somega ---- s
                %                  |       |        |
                %                  |       |        |
                %                  |       |        |
                %    SW=(i+1,j-1)  Sw - Somega  -  S=(i+1,j)
                
                
                % Calculate some help values
                
                y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
                y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
                y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
                y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
                
                y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
                y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
                y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
                
                
                
                %{
                y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
                y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
                y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
                y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
                y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
                y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
                
                y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
                
                y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
                y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                
                y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
                y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
                %}
                
                % Around somega
                
                dy_Sw_S = y_S - y_Sw;  dx_Sw_S = x_S - x_Sw;
                dy_S_P  = y_P - y_S;   dx_S_P  = x_P - x_S;
                dy_P_w  = y_w - y_P;   dx_P_w  = x_w - x_P;
                dy_w_Sw = y_Sw - y_w;  dx_w_Sw = x_Sw - x_w;
                
                % Around sigmaw
                
                dy_sW_s = y_s - y_sW;  dx_sW_s = x_s - x_sW;
                dy_s_P  = y_P - y_s;   dx_s_P  = x_P - x_s;
                dy_P_W  = y_W - y_P;   dx_P_W  = x_W - x_P;
                dy_W_sW = y_sW - y_W;  dx_W_sW = x_sW - x_W;
                
                % Around P
                
                dy_sw_s  = y_s - y_sw;  dx_sw_s  = x_s - x_sw;
                dy_w_sw  = y_sw - y_w;  dx_w_sw  = x_sw - x_w;
                
                % At boundary
                
                dl_s_P = norm([dx_s_P; dy_s_P]);
                dl_P_w = norm([dx_P_w; dy_P_w]);
                
                % Areas
                
                S_sigmaomega          = abs( x_P*y_s - y_P*x_s...
                    + x_s*y_sw - y_s*x_sw...
                    + x_sw*y_w - y_sw*x_w...
                    + x_w*y_P - y_w*x_P) / 2;
                S_somega                 = abs( x_P*y_S - y_P*x_S...
                    + x_S*y_Sw - y_S*x_Sw...
                    + x_Sw*y_w - y_Sw*x_w...
                    + x_w*y_P - y_w*x_P) / 2;
                S_sigmaw                 = abs( x_P*y_s - y_P*x_s...
                    + x_s*y_sW - y_s*x_sW...
                    + x_sW*y_W - y_sW*x_W...
                    + x_W*y_P - y_W*x_P) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_northEast
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                
                % West
                stecil(index(i, j-1)) = D_3;
                % South
                stecil(index(i+1, j)) = D1;
                % SW
                stecil(index(i+1, j-1)) = D_2;
                % P
                stecil(index(i, j)) = D0;
                
                
%             end
            %fbm......................................................
            if strcmp(str,'p')
                    if strcmp(boundary.p.east, 'neumann') || strcmp(boundary.p.east, 'wall')
                        if strcmp(boundary.p.north, 'neumann') || strcmp(boundary.p.north, 'wall')
                            b = 0;
                        elseif strcmp(boundary.p.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.p_out;
                        end
                        % please debug and correct the following section
                    elseif strcmp(boundary.p.east, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.p_out;
                    end
            elseif strcmp(str, 'u')
                    if strcmp(boundary.u.east, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.u.east, 'neumann')
                        if strcmp(boundary.u.north, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.u.north, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.u.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.u_in;
                        end
                    elseif strcmp(boundary.u.east, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.u_in;
                    end
            elseif strcmp(str, 'v')
                    if strcmp(boundary.v.east, 'wall') 
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.v.east, 'neumann')
                        if strcmp(boundary.v.north, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.v.north, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.v.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.v_in;
                        end
                    elseif strcmp(boundary.v.east, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.v_in;
                    end    
            end                  
            
            
%             if strcmp(boundary.east, 'Neumann')
%                 if strcmp(boundary.north, 'Neumann')
%                     % P
%                     stecil(index(i, j)) = D0;
%                     b=q.east*dl_s_P/S_sigmaomega + q.north*dl_P_w/S_sigmaomega;
%                 else
%                     stecil(index(i, j)) = D0 - alpha*dl_P_w/S_sigmaomega;
%                     b = -alpha*Tinf*dl_P_w/S_sigmaomega + q.east*dl_s_P/S_sigmaomega;
%                 end
%             elseif strcmp(boundary.east, 'Robin')
%                 if strcmp(boundary.north, 'Neumann')
%                     % P
%                     stecil(index(i, j)) = D0 - alpha*dl_s_P/S_sigmaomega;
%                     b = -alpha*Tinf*dl_s_P/S_sigmaomega + q.north*dl_P_w/S_sigmaomega;
%                 else
%                     stecil(index(i, j)) = D0 - (alpha*dl_P_w/S_sigmaomega)...
%                         - alpha*dl_s_P/S_sigmaomega;
%                     b = -Tinf*((alpha*dl_P_w/S_sigmaomega)...
%                         + (alpha*dl_s_P/S_sigmaomega));
%                 end
%                 end
                
                
            %fbm.......................................................    
            %Misc......................................................
%             elseif strcmp(boundary.east, 'Dirichlet')
%                 stecil(index(i, j)) = 1;
%                 b=TD.east;
%             end
%        end
            %Misc......................................................
            
            
            %        stecil(index(i,j)) = 1;
            %        if i == 1 && ~(strcmp(boundary.north, 'Dirichlet') && strcmp(boundary.east, 'Dirichlet'))
            %            stecil(index(i, j-1)) = -0.5;
            %            stecil(index(i+1, j)) = -0.5;
            %        else
            %            stecil(index(i-1, j)) = -1;
            %        end
                
%         elseif strcmp(boundary.north, 'Dirichlet')
%             stecil(index(i, j)) = 1;
%             b=TD.north;
%         end
        
        
        case 'WestSouth'
        
%         if strcmp(boundary.south, 'Neumann')
%             
%             if strcmp(boundary.west, 'Neumann') || strcmp(boundary.west, 'Robin')
                
                
                
                % Calculate some help values
                
                y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
                y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
                y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
                y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
                
                y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                
                y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2;
                y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
                
                y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
                
                % Around nepsilon
                
                dy_P_e   = y_e - y_P;    dx_P_e   = x_e - x_P;
                dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
                dy_Ne_N  = y_N - y_Ne;   dx_Ne_N  = x_N - x_Ne;
                dy_N_P   = y_P - y_N;    dx_N_P   = x_P - x_N;
                
                
                % Around etae
                
                dy_P_E  = y_E - y_P;    dx_P_E  = x_E - x_P;
                dy_E_nE = y_nE - y_E;   dx_E_nE = x_nE - x_E;
                dy_nE_n = y_n - y_nE;   dx_nE_n = x_n - x_nE;
                dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
                
                
                % Around P
                
                dy_e_ne = y_ne - y_e;  dx_e_ne = x_ne - x_e;
                dy_ne_n = y_n - y_ne;  dx_ne_n = x_n - x_ne;
                
                % At boundary
                
                dl_n_P = norm([dx_n_P; dy_n_P]);
                %fbm..............................................
                dl_P_e = norm([dx_P_e; dy_P_e]);
                %fbm..............................................
                % Areas
                
                S_etae                 = abs( x_nE*y_E - y_nE*x_E...
                    + x_E*y_P - y_E*x_P...
                    + x_P*y_n - y_P*x_n...
                    + x_n*y_nE - y_n*x_nE) / 2;
                
                
                S_nepsilon                 = abs( x_P*y_N - y_P*x_N...
                    + x_N*y_Ne - y_N*x_Ne...
                    + x_Ne*y_e - y_Ne*x_e...
                    + x_e*y_P - y_e*x_P ) / 2;
                
                
                S_etaepsilon         = abs( x_P*y_n - y_P*x_n...
                    + x_n*y_ne - y_n*x_ne...
                    + x_ne*y_e - y_ne*x_e...
                    + x_e*y_P - y_e*x_P) / 2;
                
                
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_southWest
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                % East
                stecil(index(i, j+1)) = D3;
                % North
                stecil(index(i-1, j)) = D_1;
                % NE
                stecil(index(i-1, j+1)) = D2;
%                 P
                stecil(index(i, j)) = D0;
                
                
%             end

            
            %%% NOTE: There is one face missing for the stencil in P (dl_P_e)
                if strcmp(str,'p')
                    if strcmp(boundary.p.south, 'neumann') || strcmp(boundary.p.south, 'wall')
                        if strcmp(boundary.p.west, 'neumann') || strcmp(boundary.p.west, 'wall')
                            b = 0;
                        elseif strcmp(boundary.p.west, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.p_out;
                        end
                    elseif strcmp(boundary.p.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.p_out;
                    end
                elseif strcmp(str, 'u')
                    if strcmp(boundary.u.south, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.u.south, 'neumann')
                        if strcmp(boundary.u.west, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.u.west, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.u.west, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.u_in;
                        end
                    elseif strcmp(boundary.u.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.u_in;
                    end
                elseif strcmp(str, 'v')
                    if strcmp(boundary.v.south, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.v.south, 'neumann')
                        if strcmp(boundary.v.west, 'wall')
                            stecil = zeros(1, n*m);s
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.v.west, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.v.west, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.v_in;
                        end
                    elseif strcmp(boundary.v.south, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.v_in;
                    end    
                end                
                
    case 'WestNorth'     
%         if strcmp(boundary.north, 'Neumann') || strcmp(boundary.north, 'Robin')
%             
%             if strcmp(boundary.west, 'Neumann') || strcmp(boundary.west, 'Robin')
                

                
                
                % Calculate some help values
                
                y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
                y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
                y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
                y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
                
                y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
                y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
                y_sE = (y_E  + y_SE)/2;  x_sE = (x_E  + x_SE)/2;
                y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
                y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
                
                % Around sepsilon
                
                dy_S_Se  = y_Se - y_S;   dx_S_Se = x_Se - x_S;
                dy_Se_e  = y_e - y_Se;   dx_Se_e = x_e - x_Se;
                dy_e_P   = y_P - y_e;    dx_e_P  = x_P - x_e;
                dy_P_S   = y_S - y_P;    dx_P_S  = x_S - x_P;
                
                
                % Around sigmae
                
                dy_s_sE  = y_sE - y_s ; dx_s_sE  = x_sE - x_s ;
                dy_sE_E  = y_E - y_sE ; dx_sE_E  = x_E - x_sE ;
                dy_E_P   = y_P - y_E ;  dx_E_P   = x_P - x_E ;
                dy_P_s   = y_s - y_P ;  dx_P_s   = x_s - x_P ;
                
                
                % Around P
                
                dy_s_se  = y_se - y_s;  dx_s_se  = x_se - x_s;
                dy_se_e  = y_e - y_se;  dx_se_e  = x_e - x_se;
                
                % At boundary
                
                dl_P_s = norm([dx_P_s; dy_P_s]);
                dl_e_P = norm([dx_e_P; dy_e_P]);
                
                % Areas
                
                
                S_sigmaepsilon  = abs( x_s*y_P - y_s*x_P...
                    + x_P*y_e - y_P*x_e...
                    + x_e*y_se - y_e*x_se...
                    + x_se*y_s - y_se*x_s) / 2;
                
                
                S_sepsilon                 = abs( x_S*y_P - y_S*x_P...
                    + x_P*y_e - y_P*x_e...
                    + x_e*y_Se - y_e*x_Se...
                    + x_Se*y_S - y_Se*x_S) / 2;
                
                
                S_sigmae                 = abs( x_E*y_sE - y_E*x_sE...
                    + x_sE*y_s - y_sE*x_s...
                    + x_s*y_P - y_s*x_P...
                    + x_P*y_E - y_P*x_E) / 2;
                
                
                %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$
                
                build_northWest
                
                %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                
                % East
                stecil(index(i, j+1)) = D3;
                % South
                stecil(index(i+1, j)) = D1;
                % SE
                stecil(index(i+1, j+1)) = D4;
                % P
                stecil(index(i, j)) = D0;
                
                if strcmp(str,'p')
                    if strcmp(boundary.p.west, 'neumann') || strcmp(boundary.p.west, 'wall')
                        if strcmp(boundary.p.north, 'neumann') || strcmp(boundary.p.north, 'wall')
                            b = 0;
                        elseif strcmp(boundary.p.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.p_out;
                        end
                    elseif strcmp(boundary.p.west, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = value.p_out;
                    end
                elseif strcmp(str, 'u')
                    if strcmp(boundary.u.west, 'wall')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.u.west, 'neumann')
                        if strcmp(boundary.u.north, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.u.north, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.u.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        end
                    elseif strcmp(boundary.u.west, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    end
                elseif strcmp(str, 'v')
                    if strcmp(boundary.v.west, 'wall') 
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    elseif strcmp(boundary.v.west, 'neumann')
                        if strcmp(boundary.v.north, 'wall')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = 0;
                        elseif strcmp(boundary.v.north, 'neumann')
                            b = 0;
                        elseif strcmp(boundary.v.north, 'dirichlet')
                            stecil = zeros(1, n*m);
                            stecil(index(i, j)) = 1;
                            b = value.v_in;
                        end
                    elseif strcmp(boundary.v.west, 'dirichlet')
                        stecil = zeros(1, n*m);
                        stecil(index(i, j)) = 1;
                        b = 0;
                    end
                end      

        
end