function [Dx,Dy,Area,Length] = buildDxDy(space)
% This function builds the Dx and Dy matrix operator is based on the finite
% volume formulation of Camilo (see lecture script and videos).

% Initialize
Dx = zeros(space.dimY*space.dimX);
Dy= zeros(space.dimY*space.dimX);
Area.S_P = zeros(space.dimY*space.dimX,1);
Area.S_s = zeros(space.dimY*space.dimX,1);
Area.S_e = zeros(space.dimY*space.dimX,1);
Area.S_n = zeros(space.dimY*space.dimX,1);
Area.S_w = zeros(space.dimY*space.dimX,1);

Length.dy_sw_se = zeros(space.dimY*space.dimX,1);   
Length.dx_sw_se = zeros(space.dimY*space.dimX,1);
Length.dy_se_ne = zeros(space.dimY*space.dimX,1);   
Length.dx_se_ne = zeros(space.dimY*space.dimX,1);
Length.dy_ne_nw = zeros(space.dimY*space.dimX,1);   
Length.dx_ne_nw = zeros(space.dimY*space.dimX,1);
Length.dy_nw_sw = zeros(space.dimY*space.dimX,1);   
Length.dx_nw_sw = zeros(space.dimY*space.dimX,1);
%Index help function
index = @(ii, jj) ii + (jj-1) * space.dimY;
% Fill matrix for inner nodes
for i = 2:space.dimY-1
    for j = 2:space.dimX-1
        % Get all the necessary geometric properties
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
        y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
        y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
        y_nE =  (y_NE + y_E)/2;  x_nE =  (x_NE + x_E)/2;
        y_sE =  (y_E +  y_SE)/2; x_sE =  (x_E +  x_SE)/2;
        y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;
        y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
        y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
        y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
        % Lengths around P
        dy_sw_se = y_se - y_sw;   dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;   dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;   dx_nw_sw = x_sw - x_nw;
        % Area
        S_P                 = abs( x_ne*y_se - y_ne*x_se...
            + x_se*y_sw - y_se*x_sw...
            + x_sw*y_nw - y_sw*x_nw...
            + x_nw*y_ne - y_nw*x_ne) / 2;
        Length.dy_sw_se(index(i,j)) = y_se - y_sw;   
        Length.dx_sw_se(index(i,j)) = x_se - x_sw;
        Length.dy_se_ne(index(i,j)) = y_ne - y_se;   
        Length.dx_se_ne(index(i,j)) = x_ne - x_se;
        Length.dy_ne_nw(index(i,j)) = y_nw - y_ne;   
        Length.dx_ne_nw(index(i,j)) = x_nw - x_ne;
        Length.dy_nw_sw(index(i,j)) = y_sw - y_nw;   
        Length.dx_nw_sw(index(i,j)) = x_sw - x_nw;
        % Areas
        Area.S_P(index(i,j))                 = abs( x_ne*y_se - y_ne*x_se...
            + x_se*y_sw - y_se*x_sw...
            + x_sw*y_nw - y_sw*x_nw...
            + x_nw*y_ne - y_nw*x_ne) / 2;
        Area.S_s(index(i+1,j))                 = abs( x_e*y_Se - y_e*x_Se...
            + x_Se*y_Sw - y_Se*x_Sw...
            + x_Sw*y_w - y_Sw*x_w...
            + x_w*y_e - y_w*x_e) / 2;
        Area.S_e(index(i,j+1))                 = abs( x_nE*y_sE - y_nE*x_sE...
            + x_sE*y_s - y_sE*x_s...
            + x_s*y_n - y_s*x_n...
            + x_n*y_nE - y_n*x_nE) / 2;
        Area.S_n(index(i-1,j))                 = abs( x_Ne*y_e - y_Ne*x_e...
            + x_e*y_w - y_e*x_w...
            + x_w*y_Nw - y_w*x_Nw...
            + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
        Area.S_w(index(i,j-1))                 = abs( x_n*y_s - y_n*x_s...
            + x_s*y_sW - y_s*x_sW...
            + x_sW*y_nW - y_sW*x_nW...
            + x_nW*y_n - y_nW*x_n) / 2;
        %Build matrix
        %South
        Dx(index(i,j),index(i+1,j)) = 0.5*dy_sw_se/S_P;
        Dy(index(i,j),index(i+1,j)) = -0.5*dx_sw_se/S_P;
        %East
        Dx(index(i,j),index(i,j+1)) = 0.5*dy_se_ne/S_P;
        Dy(index(i,j),index(i,j+1)) = -0.5*dx_se_ne/S_P;
        %North
        Dx(index(i,j),index(i-1,j)) = 0.5*dy_ne_nw/S_P;
        Dy(index(i,j),index(i-1,j)) = -0.5*dx_ne_nw/S_P;
        %West
        Dx(index(i,j),index(i,j-1)) = 0.5*dy_nw_sw/S_P;
        Dy(index(i,j),index(i,j-1)) = -0.5*dx_nw_sw/S_P;
        %Node itself
        Dx(index(i,j),index(i,j)) = (0.5*dy_sw_se + 0.5*dy_se_ne + 0.5*dy_ne_nw + 0.5*dy_nw_sw)/S_P;
        Dy(index(i,j),index(i,j)) = -(0.5*dx_sw_se + 0.5*dx_se_ne + 0.5*dx_ne_nw + 0.5*dx_nw_sw)/S_P;
    end
end

% Note: Currently, Dx and Dy are only calculated for the inner nodes.
% However, the commented lines below allow to calculate the boundary values
% based on Finite Volumes using only information from the inner nodes and
% artificial, help integration points. If no information on Dx and Dy
% is available at the boundaries, you can potentially use the following 
% lines.

% % South boundary
% i=space.dimY;
% for j=2:space.dimX-1
%     y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
%     y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
%     y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
%     y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
%     y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
%     y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
%     y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
%     y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
%     y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
%     y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
%     y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
%     y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;
%     % Around n
%     dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;  
%     % Around P
%     dy_e_ne = y_ne - y_e;     dx_e_ne = x_ne - x_e;
%     dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
%     dy_nw_w =  y_w - y_nw;    dx_nw_w =  x_w - x_nw;
%     % Area
%     S_eta                 = abs( x_ne*y_e - y_ne*x_e...
%         + x_e*y_w - y_e*x_w...
%         + x_w*y_nw - y_w*x_nw...
%         + x_nw*y_ne - y_nw*x_ne) / 2;
%     %East Node
%     Dx(index(i,j),index(i,j+1)) = ((3/8)*dy_e_ne)/S_eta;
%     Dy(index(i,j),index(i,j+1)) = -((3/8)*dx_e_ne)/S_eta;
%     %North Node
%     Dx(index(i,j),index(i-1,j)) = ((1/8)*dy_e_ne + 0.5*dy_ne_nw + (1/8)*dy_nw_w)/S_eta;
%     Dy(index(i,j),index(i-1,j)) = -((1/8)*dx_e_ne + 0.5*dx_ne_nw + (1/8)*dx_nw_w)/S_eta;
%     %West Node
%     Dx(index(i,j),index(i,j-1)) = ((3/8)*dy_nw_w)/S_eta;
%     Dy(index(i,j),index(i,j-1)) = -((3/8)*dx_nw_w)/S_eta;
%     %Node itself
%     Dx(index(i,j),index(i,j)) = (dy_w_e + (3/8)*dy_e_ne + 0.5*dy_ne_nw + (3/8)*dy_nw_w)/S_eta;
%     Dy(index(i,j),index(i,j)) = -(dx_w_e + (3/8)*dx_e_ne + 0.5*dx_ne_nw + (3/8)*dx_nw_w)/S_eta;
%     %Northeast Node
%     Dx(index(i,j),index(i-1,j+1)) = ((1/8)*dy_e_ne)/S_eta;
%     Dy(index(i,j),index(i-1,j+1)) = -((1/8)*dx_e_ne)/S_eta;
%     %Northwest Node
%     Dx(index(i,j),index(i-1,j-1)) = ((1/8)*dy_nw_w)/S_eta;
%     Dy(index(i,j),index(i-1,j-1)) = -((1/8)*dx_nw_w)/S_eta;
% end
% 
% % North boundary
% i = 1;
% for j=2:space.dimX-1
%     y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
%     y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
%     y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
%     y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
%     y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
%     y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
%     y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
%     y_e  = (y_E  + y_P)/2;  x_e  = (x_E  + x_P)/2;
%     y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
%     y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
%     y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
%     y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
%     % Around P
%     dy_sw_se  = y_se - y_sw;  dx_sw_se  = x_se - x_sw;
%     dy_se_e   = y_e - y_se ;  dx_se_e   = x_e - x_se ;
%     dy_e_w    = y_w - y_e;    dx_e_w    = x_w - x_e;
%     dy_w_sw   = y_sw - y_w;   dx_w_sw   = x_sw - x_w;
%     % Areas
%     S_sigma                 = abs( x_e*y_se - y_e*x_se...
%         + x_se*y_sw - y_se*x_sw...
%         + x_sw*y_w - y_sw*x_w...
%         + x_w*y_e - y_w*x_e) / 2;
%     %East Node
%     Dx(index(i,j),index(i,j+1)) = ((3/8)*dy_se_e)/S_sigma;
%     Dy(index(i,j),index(i,j+1)) = -((3/8)*dx_se_e)/S_sigma;
%     %South Node
%     Dx(index(i,j),index(i+1,j)) = ((1/8)*dy_se_e + 0.5*dy_sw_se + (1/8)*dy_w_sw)/S_sigma;
%     Dy(index(i,j),index(i+1,j)) = -((1/8)*dx_se_e + 0.5*dx_sw_se + (1/8)*dx_w_sw)/S_sigma;
%     %West Node
%     Dx(index(i,j),index(i,j-1)) = ((3/8)*dy_w_sw)/S_sigma;
%     Dy(index(i,j),index(i,j-1)) = -((3/8)*dx_w_sw)/S_sigma;
%     %Node itself
%     Dx(index(i,j),index(i,j)) = (dy_e_w + (3/8)*dy_se_e + 0.5*dy_sw_se + (3/8)*dy_w_sw)/S_sigma;
%     Dy(index(i,j),index(i,j)) = -(dx_e_w + (3/8)*dx_se_e + 0.5*dx_sw_se + (3/8)*dx_w_sw)/S_sigma;
%     %Southeast Node
%     Dx(index(i,j),index(i+1,j+1)) = ((1/8)*dy_se_e)/S_sigma;
%     Dy(index(i,j),index(i+1,j+1)) = -((1/8)*dx_se_e)/S_sigma;
%     %Southwest Node
%     Dx(index(i,j),index(i+1,j-1)) = ((1/8)*dy_w_sw)/S_sigma;
%     Dy(index(i,j),index(i+1,j-1)) = -((1/8)*dx_w_sw)/S_sigma;
% end
% 
% % West boundary
% j=1;
% for i=2:space.dimY-1
%     y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
%     y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
%     y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
%     y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
%     y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
%     y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  ); 
%     y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
%     y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
%     y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2; 
%     y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
%     y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
%     y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
%     y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
%     % Around P
%     dy_s_se  = y_se - y_s;   dx_s_se  = x_se - x_s;
%     dy_se_ne = y_ne - y_se;  dx_se_ne = x_ne - x_se;
%     dy_ne_n  = y_n - y_ne;   dx_ne_n  = x_n - x_ne;
%     dy_n_s   = y_s - y_n;    dx_n_s   = x_s - x_n;
%     % Area
%     S_epsilon                 = abs( x_s*y_n - y_s*x_n...
%         + x_n*y_ne - y_n*x_ne...
%         + x_ne*y_se - y_ne*x_se...
%         + x_se*y_s - y_se*x_s) / 2;
%     %South
%     Dx(index(i,j),index(i+1,j)) = ((3/8)*dy_s_se)/S_epsilon;
%     Dy(index(i,j),index(i+1,j)) = -((3/8)*dx_s_se)/S_epsilon;
%     %East
%     Dx(index(i,j),index(i,j+1)) = ((1/8)*dy_s_se + 0.5*dy_se_ne + (1/8)*dy_ne_n)/S_epsilon;
%     Dy(index(i,j),index(i,j+1)) = -((1/8)*dx_s_se + 0.5*dx_se_ne + (1/8)*dx_ne_n)/S_epsilon;
%     %North
%     Dx(index(i,j),index(i-1,j)) = ((3/8)*dy_ne_n)/S_epsilon;
%     Dy(index(i,j),index(i-1,j)) = -((3/8)*dx_ne_n)/S_epsilon;
%     %Node itself
%     Dx(index(i,j),index(i,j)) = (dy_n_s + (3/8)*dy_s_se + (3/8)*dy_ne_n + 0.5*dy_se_ne)/S_epsilon;
%     Dy(index(i,j),index(i,j)) = -(dx_n_s + (3/8)*dx_s_se + (3/8)*dx_ne_n + 0.5*dx_se_ne)/S_epsilon;
%     %Northeast
%     Dx(index(i,j),index(i-1,j+1)) = (1/8)*dy_ne_n/S_epsilon;
%     Dy(index(i,j),index(i-1,j+1)) = -(1/8)*dx_ne_n/S_epsilon;
%     %Southeast
%     Dx(index(i,j),index(i+1,j+1)) = (1/8)*dy_s_se/S_epsilon;
%     Dy(index(i,j),index(i+1,j+1)) = -(1/8)*dx_s_se/S_epsilon;
% end
% 
% % East boundary
% j=space.dimX;
% for i=2:space.dimY-1
%     y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
%     y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
%     y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
%     y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
%     y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
%     y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
%     y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
%     y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
%     y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
%     y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
%     y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
%     y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
%     y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
%     % Around P
%     dy_sw_s  = y_s - y_sw;   dx_sw_s  = x_s - x_sw;
%     dy_s_n   = y_n - y_s;    dx_s_n   = x_n - x_s;
%     dy_n_nw  = y_nw - y_n;   dx_n_nw  = x_nw - x_n;
%     dy_nw_sw = y_sw - y_nw;  dx_nw_sw = x_sw - x_nw;
%     % Area
%     S_omega                 = abs( x_n*y_s - y_n*x_s...
%         + x_s*y_sw - y_s*x_sw...
%         + x_sw*y_nw - y_sw*x_nw...
%         + x_nw*y_n - y_nw*x_n) / 2;
%     %South
%     Dx(index(i,j),index(i+1,j)) = ((3/8)*dy_sw_s)/S_omega;
%     Dy(index(i,j),index(i+1,j)) = -((3/8)*dx_sw_s)/S_omega;
%     %West
%     Dx(index(i,j),index(i,j-1)) = ((1/8)*dy_sw_s + 0.5*dy_nw_sw + (1/8)*dy_n_nw)/S_omega;
%     Dy(index(i,j),index(i,j-1)) = -((1/8)*dx_sw_s + 0.5*dx_nw_sw + (1/8)*dx_n_nw)/S_omega;
%     %North
%     Dx(index(i,j),index(i-1,j)) = ((3/8)*dy_n_nw)/S_omega;
%     Dy(index(i,j),index(i-1,j)) = -((3/8)*dx_n_nw)/S_omega;
%     %Node itself
%     Dx(index(i,j),index(i,j)) = (dy_s_n + (3/8)*dy_sw_s + (3/8)*dy_n_nw + 0.5*dy_nw_sw)/S_omega;
%     Dy(index(i,j),index(i,j)) = -(dx_s_n + (3/8)*dx_sw_s + (3/8)*dx_n_nw + 0.5*dx_nw_sw)/S_omega;
%     %Northwest
%     Dx(index(i,j),index(i-1,j-1)) = (1/8)*dy_n_nw/S_omega;
%     Dy(index(i,j),index(i-1,j-1)) = -(1/8)*dx_n_nw/S_omega;
%     %Southwest
%     Dx(index(i,j),index(i+1,j-1)) = (1/8)*dy_sw_s/S_omega;
%     Dy(index(i,j),index(i+1,j-1)) = -(1/8)*dx_sw_s/S_omega;
% end
% 
% % Southeast Corner
% i = space.dimY;
% j = space.dimX;
% y_NW = space.Y(i-1,j-1);   x_NW = space.X(i-1,j-1);
% y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
% y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
% y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
% y_Nw = (y_N + y_NW)/2;  x_Nw = (x_N + x_NW)/2;
% y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
% y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
% y_nw = (y_w + y_Nw)/2;  x_nw = (x_w + x_Nw)/2;
% % Around nomega
% dy_w_P  = y_P - y_w;   dx_w_P  = x_P - x_w;
% % Around etaw
% dy_P_n  = y_n - y_P;    dx_P_n  = x_n - x_P;
% % Around P
% dy_n_nw = y_nw - y_n;  dx_n_nw = x_nw - x_n;
% dy_nw_w = y_w - y_nw;  dx_nw_w = x_w - x_nw;
% % Area
% S_etaomega         = abs( x_n*y_P - y_n*x_P...
%     + x_P*y_w - y_P*x_w...
%     + x_w*y_nw - y_w*x_nw...
%     + x_nw*y_n - y_nw*x_n) / 2;
% %North node
% Dx(index(i,j),index(i-1,j)) = (0.25*dy_P_n+(3/8)*dy_n_nw+(1/8)*dy_nw_w)/S_etaomega;
% Dy(index(i,j),index(i-1,j)) = -(0.25*dx_P_n+(3/8)*dx_n_nw+(1/8)*dx_nw_w)/S_etaomega;
% %West node
% Dx(index(i,j),index(i,j-1)) = ((1/8)*dy_n_nw+(3/8)*dy_nw_w+0.25*dy_w_P)/S_etaomega;
% Dy(index(i,j),index(i,j-1)) = -((1/8)*dx_n_nw+(3/8)*dx_nw_w+0.25*dx_w_P)/S_etaomega;
% %Northwest node
% Dx(index(i,j),index(i-1,j-1)) = ((1/8)*dy_n_nw+(1/8)*dy_nw_w)/S_etaomega;
% Dy(index(i,j),index(i-1,j-1)) = -((1/8)*dx_n_nw+(1/8)*dx_nw_w)/S_etaomega;
% %Node itself
% Dx(index(i,j),index(i,j)) = (0.75*dy_P_n+(3/8)*dy_n_nw+(3/8)*dy_nw_w+0.75*dy_w_P)/S_etaomega;
% Dy(index(i,j),index(i,j)) = -(0.75*dx_P_n+(3/8)*dx_n_nw+(3/8)*dx_nw_w+0.75*dx_w_P)/S_etaomega;
% 
% % Northeast corner
% i=1;
% j=space.dimX;
% y_W  = space.Y(i  ,j-1);   x_W  = space.X(i  ,j-1);
% y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
% y_SW = space.Y(i+1,j-1);   x_SW = space.X(i+1,j-1);
% y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
% y_w  = (y_P  + y_W)/2;  x_w  = (x_P  + x_W)/2;
% y_Sw = (y_S  + y_SW)/2; x_Sw = (x_S  + x_SW)/2;
% y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
% y_sw = (y_w + y_Sw)/2;  x_sw = (x_w + x_Sw)/2;
% % Around somega
% dy_P_w  = y_w - y_P;   dx_P_w  = x_w - x_P;
% % Around sigmaw
% dy_s_P  = y_P - y_s;   dx_s_P  = x_P - x_s;
% % Around P
% dy_sw_s  = y_s - y_sw;  dx_sw_s  = x_s - x_sw;
% dy_w_sw  = y_sw - y_w;  dx_w_sw  = x_sw - x_w;
% % Area
% S_sigmaomega          = abs( x_P*y_s - y_P*x_s...
%     + x_s*y_sw - y_s*x_sw...
%     + x_sw*y_w - y_sw*x_w...
%     + x_w*y_P - y_w*x_P) / 2;
% %West node
% Dx(index(i,j),index(i,j-1)) = (0.25*dy_P_w+(3/8)*dy_w_sw+(1/8)*dy_sw_s)/S_sigmaomega;
% Dy(index(i,j),index(i,j-1)) = -(0.25*dx_P_w+(3/8)*dx_w_sw+(1/8)*dx_sw_s)/S_sigmaomega;
% %South node
% Dx(index(i,j),index(i+1,j)) = ((1/8)*dy_w_sw+(3/8)*dy_sw_s+0.25*dy_s_P)/S_sigmaomega;
% Dy(index(i,j),index(i+1,j)) = -((1/8)*dx_w_sw+(3/8)*dx_sw_s+0.25*dx_s_P)/S_sigmaomega;
% %Southwest node
% Dx(index(i,j),index(i+1,j-1)) = ((1/8)*dy_w_sw+(1/8)*dy_sw_s)/S_sigmaomega;
% Dy(index(i,j),index(i+1,j-1)) = -((1/8)*dx_w_sw+(1/8)*dx_sw_s)/S_sigmaomega;
% %Node itself
% Dx(index(i,j),index(i,j)) = (0.75*dy_P_w+(3/8)*dy_w_sw+(3/8)*dy_sw_s+0.75*dy_s_P)/S_sigmaomega;
% Dy(index(i,j),index(i,j)) = -(0.75*dx_P_w+(3/8)*dx_w_sw+(3/8)*dx_sw_s+0.75*dx_s_P)/S_sigmaomega;
% 
% % Northwest corner
% i=1;
% j=1;
% y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
% y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
% y_SE = space.Y(i+1,j+1);   x_SE = space.X(i+1,j+1);
% y_S  = space.Y(i+1,j  );   x_S  = space.X(i+1,j  );
% y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
% y_Se = (y_S  + y_SE)/2; x_Se = (x_S  + x_SE)/2;
% y_s  =  (y_P +  y_S)/2;  x_s  =  (x_P +  x_S)/2;
% y_se = (y_e + y_Se)/2;  x_se = (x_e + x_Se)/2;
% % Around sepsilon
% dy_e_P   = y_P - y_e;    dx_e_P  = x_P - x_e;
% % Around sigmae
% dy_P_s   = y_s - y_P ;  dx_P_s   = x_s - x_P ;
% % Around P
% dy_s_se  = y_se - y_s;  dx_s_se  = x_se - x_s;
% dy_se_e  = y_e - y_se;  dx_se_e  = x_e - x_se;
% % Area
% S_sigmaepsilon  = abs( x_s*y_P - y_s*x_P...
%     + x_P*y_e - y_P*x_e...
%     + x_e*y_se - y_e*x_se...
%     + x_se*y_s - y_se*x_s) / 2;
% %South node
% Dx(index(i,j),index(i+1,j)) = (0.25*dy_P_s+(3/8)*dy_s_se+(1/8)*dy_se_e)/S_sigmaepsilon;
% Dy(index(i,j),index(i+1,j)) = -(0.25*dx_P_s+(3/8)*dx_s_se+(1/8)*dx_se_e)/S_sigmaepsilon;
% %East node
% Dx(index(i,j),index(i,j+1)) = ((1/8)*dy_s_se+(3/8)*dy_se_e+0.25*dy_e_P)/S_sigmaepsilon;
% Dy(index(i,j),index(i,j+1)) = -((1/8)*dx_s_se+(3/8)*dx_se_e+0.25*dx_e_P)/S_sigmaepsilon;
% %Southeast node
% Dx(index(i,j),index(i+1,j+1)) = ((1/8)*dy_s_se+(1/8)*dy_se_e)/S_sigmaepsilon;
% Dy(index(i,j),index(i+1,j+1)) = -((1/8)*dx_s_se+(1/8)*dx_se_e)/S_sigmaepsilon;
% %Node itself
% Dx(index(i,j),index(i,j)) = (0.75*dy_P_s+(3/8)*dy_s_se+(3/8)*dy_se_e+0.75*dy_e_P)/S_sigmaepsilon;
% Dy(index(i,j),index(i,j)) = -(0.75*dx_P_s+(3/8)*dx_s_se+(3/8)*dx_se_e+0.75*dx_e_P)/S_sigmaepsilon;
% 
% % Southwest corner
% i=space.dimY;
% j=1;
% y_NE = space.Y(i-1,j+1);   x_NE = space.X(i-1,j+1);
% y_N  = space.Y(i-1,j  );   x_N  = space.X(i-1,j  );
% y_E  = space.Y(i  ,j+1);   x_E  = space.X(i  ,j+1);
% y_P  = space.Y(i  ,j  );   x_P  = space.X(i  ,j  );
% y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
% y_e  = (y_P  + y_E)/2;  x_e  = (x_P  + x_E)/2;
% y_n  =  (y_N + y_P)/2;   x_n  =  (x_N + x_P)/2;
% y_ne = (y_e + y_Ne)/2;  x_ne = (x_e + x_Ne)/2;
% % Around nepsilon
% dy_P_e   = y_e - y_P;    dx_P_e   = x_e - x_P;
% % Around etae
% dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
% % Around P
% dy_e_ne = y_ne - y_e;  dx_e_ne = x_ne - x_e;
% dy_ne_n = y_n - y_ne;  dx_ne_n = x_n - x_ne;
% % Area
% S_etaepsilon         = abs( x_P*y_n - y_P*x_n...
%     + x_n*y_ne - y_n*x_ne...
%     + x_ne*y_e - y_ne*x_e...
%     + x_e*y_P - y_e*x_P) / 2;
% %North node
% Dx(index(i,j),index(i-1,j)) = ((1/8)*dy_e_ne+(3/8)*dy_ne_n+0.25*dy_n_P)/S_etaepsilon;
% Dy(index(i,j),index(i-1,j)) = -((1/8)*dx_e_ne+(3/8)*dx_ne_n+0.25*dx_n_P)/S_etaepsilon;
% %East node
% Dx(index(i,j),index(i,j+1)) = (0.25*dy_P_e+(3/8)*dy_e_ne+(1/8)*dy_ne_n)/S_etaepsilon;
% Dy(index(i,j),index(i,j+1)) = -(0.25*dx_P_e+(3/8)*dx_e_ne+(1/8)*dx_ne_n)/S_etaepsilon;
% %Northeast node
% Dx(index(i,j),index(i-1,j+1)) = ((1/8)*dy_e_ne+(1/8)*dy_ne_n)/S_etaepsilon;
% Dy(index(i,j),index(i-1,j+1)) = -((1/8)*dx_e_ne+(1/8)*dx_ne_n)/S_etaepsilon;
% %Node itself
% Dx(index(i,j),index(i,j)) = (0.75*dy_P_e+(3/8)*dy_e_ne+(3/8)*dy_ne_n+0.75*dy_n_P)/S_etaepsilon;
% Dy(index(i,j),index(i,j)) = -(0.75*dx_P_e+(3/8)*dx_e_ne+(3/8)*dx_ne_n+0.75*dx_n_P)/S_etaepsilon;

% Make results sparse to save space
Dx = sparse(Dx);
Dy = sparse(Dy);

end