%% Fluid Variables
FV.rho = 1.177;               % Dry air at 300K
FV.mu = 1.846e-02;            % Dynamic viscosity (altered for Re)
FV.nu = FV.mu/FV.rho;         % Kinematic viscosity

%% Time Variables and Numerical Scheme
time.tol = 1e-2;              % Tolerance criterium for continuity eq.
time.iter = 200;              % Purpose unknown yet
time.dt = 0.01;               % Time step size

%% Domain Discretization and Mesh

space.l = 1;                  % Channel length
space.h = 0.2;                % Channel height
space.dimX = 250;             % Number of nodes in x direction
space.dimY = 50;              % Number of nodes in y direction

% Generate Mesh
x = linspace(0,space.l,space.dimX);
space.X = meshgrid(x,linspace(0,space.h,space.dimY));
space.Y = zeros(space.dimY,space.dimX);

% Mesh generation for plain channel flow
for i=1:1:space.dimX
   space.Y(:,i) = linspace(space.h,0,space.dimY);
end

% Formfunctions for obstacle generation

%1) Exponential obstacle

% mag = 0.4*space.h;
% pos = 0.4;
% epsilon = 50;
% formfunction = @(xnorm) mag*exp((-epsilon*(xnorm-pos)^2));
% for i=1:1:space.dimX
%    space.Y(:,i) = linspace(space.h,formfunction(x(i)),space.dimY);
% end
 
%2) Sinus obstacle
% rad = 0.5*space.h;
% pos = 0.3;

% scale = 1;
% [val,idx] = min(abs(x-(pos-rad)));
% n = idx;
% [val,idx] = min(abs(x-(pos+rad)));
% m = idx;
% theta = [zeros(1,n-1) linspace(0,pi,length(x(n:m))) zeros(1,space.dimX-m+1)];
% formfunction = @(xnorm) scale*rad*sin(xnorm);
% for i=1:1:space.dimX
%    space.Y(:,i) = linspace(space.h,formfunction(theta(i)),space.dimY);
% end

%3) Parabolic obstacle

% h = 0.5;
% k = 0.1;
% a = -40;
% e = 0.04;
% [val,idx] = min(abs(x-(h-2*e)));
% n = idx;
% [val,idx] = min(abs(x-(h+2*e)));
% m = idx; 
% t = [zeros(1,n-1) linspace(-e,e,length(x(n:m))) zeros(1,space.dimX-m+1)];
% % t = linspace(-e,e,length(x))
% % formfunction = @(xnorm) a*t.^2 + k;
% for i=1:1:space.dimX
%    space.Y(:,i) = linspace(space.h,form(a,t(i),k,h),space.dimY);
% end

%4) Circular obstacle
% pos = 0.4;
% r = 0.1;
% formfunction = @(xnorm) real(sqrt(r^2-(xnorm-pos)^2));
% for i=1:1:space.dimX
%    space.Y(:,i) = linspace(space.h,formfunction(x(i)),space.dimY);
% end

% Initialize Matrices for pressure and velocity
u = zeros(space.dimY,space.dimX);
v = zeros(space.dimY,space.dimX);
p = zeros(space.dimY,space.dimX);

%% Boundary conditions
% In order to solve the Navier-Stokes equations, boundary conditions for
% both p and u are necessary. For this test case, we model a channel flow
% from left to right which requires velocity definition at the inlet and
% pressure definition at the outlet.

boundary.u.north = 'wall';
boundary.u.south = 'wall';
boundary.u.west = 'dirichlet';
boundary.u.east = 'neumann';

boundary.v.north = 'wall';
boundary.v.south = 'wall';
boundary.v.west = 'dirichlet';
boundary.v.east = 'neumann';

boundary.p.north = 'neumann';
boundary.p.south = 'neumann';
boundary.p.west = 'neumann';   
boundary.p.east = 'dirichlet'; 

% Values for the inlet and outlet must be specified
value.u_in = 0.1;       % u inlet velocity 0.1 m/s
value.v_in = 0;         % v inlet velocity 0 m/s
value.p_out = 0;        % Outlet pressure 0 Pa

% Initialize logical arrays
ib.u = true(space.dimY,space.dimX);
ib.v = true(space.dimY,space.dimX);
ib.pre = true(space.dimY,space.dimX);
ib.rhs = true(space.dimY,space.dimX);

ib.P = false(space.dimY,space.dimX);
ib.s = false(space.dimY,space.dimX);
ib.e = false(space.dimY,space.dimX);
ib.n = false(space.dimY,space.dimX);
ib.w = false(space.dimY,space.dimX);

ib.P(2:space.dimY-1,2:space.dimX-1) = true;
ib.s(3:space.dimY,2:space.dimX-1) = true;
ib.e(2:space.dimY-1,3:space.dimX) = true;
ib.n(1:space.dimY-2,2:space.dimX-1) = true;
ib.w(2:space.dimY-1,1:space.dimX-2) = true;

ib.P = ib.P(:);
ib.s = ib.s(:);
ib.e = ib.e(:);
ib.n = ib.n(:);
ib.w = ib.w(:);

ib.u(1,1:space.dimX) = false;
ib.u(2:space.dimY-1,1) = false;
ib.u(space.dimY,1:space.dimX) = false;
ib.v(1,1:space.dimX) = false;
ib.v(1:space.dimY,1) = false;
ib.v(space.dimY,1:space.dimX) = false;
ib.v(1:space.dimY,space.dimX) = false;
ib.pre(1:space.dimY, space.dimX) = false;
ib.rhs(1,1:space.dimX) = false;
ib.rhs(1:space.dimY,1) = false;
ib.rhs(space.dimY,1:space.dimX) = false;
ib.rhs(1:space.dimY,space.dimX) = false;

% Calculate and print CFL number
CFL = (value.u_in*time.dt)/(space.l/(space.dimX-1));
message = ['The CFL number is currently ',num2str(CFL)];
disp(message);
disp('In case convergence issues occur, please reduce CFL');

% Help function for parabolic obstacle generation
function [y] = form(a,xnorm,k,h)
    if xnorm == 0.0
        y = 0;
    else
%         x = 2*xnorm*a + h;
        y = (a*xnorm^2) + k;
%         y = ((x - h)^2)/4*a + k;
    end
end

%% Initial guess for velocity and pressure
% In case the user wants to make a guess for u and p in the domain, the
% corresponding matrices can be manipulated in the section below

% No initial guesses for now

%% End of problem initialization