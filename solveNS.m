function [u,v,p,res] = solveNS(u,v,p,Dx,Dy,Lu,Lv,Lp,Bu,Bv,Bp,FV,time,space,value,ib,Area,Length)

% Initialize vector to store residuals
res = zeros(1,time.iter);

% Index help function
index = @(ii, jj) ii + (jj-1) * size(space.X, 1);

% Reshape matrices to vector
u_vec = reshape(u,[space.dimY*space.dimX,1]);
v_vec = reshape(v,[space.dimY*space.dimX,1]);
p_vec = reshape(p,[space.dimY*space.dimX,1]);

% Flatten logical arrays
ib.u = ib.u(:);
ib.v = ib.v(:);
BD = [ib.u;ib.v];
ib.pre = ib.pre(:);
ib.rhs = ib.rhs(:);

% Initialize matrices and set boundary conditions
pc = zeros(space.dimY*space.dimX,1);
pnew = zeros(space.dimY*space.dimX,1);
unew = zeros(space.dimY*space.dimX,1);
vnew = zeros(space.dimY*space.dimX,1);
u_vec(~ib.u) = Bu(~ib.u);

unew(~ib.u) = Bu(~ib.u);
v_vec(~ib.v) = Bv(~ib.v);
vnew(~ib.v) = Bv(~ib.v);
p_vec(~ib.pre) = Bp(~ib.pre);
pnew(~ib.pre) = Bp(~ib.pre);
star = [u_vec;v_vec];
pc = p_vec;
rhs = zeros(space.dimY*space.dimX,1);

% Main time loop
for t=1:time.iter
    
    % 1. Compute u* and p*
    
    % Build matrices
    A_U_U = (1/time.dt) * speye(space.dimY*space.dimX,space.dimY*space.dimX)...
        + spdiags(Dx*u_vec,0,space.dimY*space.dimX,space.dimY*space.dimX) - FV.nu*Lu;
    A_V_U = spdiags(Dy*u_vec,0,space.dimY*space.dimX,space.dimY*space.dimX);
    b_u = u_vec/time.dt - (1/FV.rho)*Dx*p_vec;
    A_U_V = spdiags(Dx*v_vec,0,space.dimY*space.dimX,space.dimY*space.dimX);
    A_V_V = (1/time.dt) * speye(space.dimY*space.dimX,space.dimY*space.dimX)...
        + spdiags(Dx*v_vec,0,space.dimY*space.dimX,space.dimY*space.dimX)-FV.nu*Lv;
    b_v = v_vec/time.dt - (1/FV.rho)*Dy*p_vec;
    
    % Assemble complete linear system
    A = [A_U_U A_V_U;A_U_V A_V_V];
    b = [b_u;b_v];
    
    % Solve system and extract u* and v*
    A = sparse(A);
    star = A\b;
    ustar = star(1:space.dimY*space.dimX);
    vstar = star(space.dimY*space.dimX+1:end);
    
    % 2. Rhie-Chow interpolation for pressure-correction equation
    [Dxustar,Dyvstar] = DxDyRC(A_U_U,A_V_V,Dx,Dy,ustar,vstar,p_vec,space,FV,ib,Area,Length);
    
    % If you do not want RC interpolation, comment out the line above and
    % uncomment the two following commands
    
    %Dxustar = Dx*ustar;
    %Dyvstar = Dy*vstar;
    
    % 3. Apply pressure correction to get p_c
    rhs(ib.rhs) = (FV.rho/time.dt)*(Dxustar(ib.rhs) + Dyvstar(ib.rhs));
    rhs(~ib.rhs) = Bp(~ib.rhs);
    
    % Solve linear system for pressure pc
    pc(ib.pre) = Lp(ib.pre,ib.pre)\rhs(ib.pre);
    
    % 4. Determine velocities u_n+1 and v_n+1 
    Dxpc = Dx*pc;
    Dypc = Dy*pc;
    
    unew(ib.u) = ustar(ib.u) - (time.dt/FV.rho)*(Dxpc(ib.u));
    vnew(ib.v) = vstar(ib.v) - (time.dt/FV.rho)*(Dypc(ib.v));
    
    u_vec(ib.u) = unew(ib.u);
    v_vec(ib.v) = vnew(ib.v);

    % 5. Update pressure p_n+1
    pnew(ib.pre) = p_vec(ib.pre) + pc(ib.pre);
    
    % Set pressure as starting point for next iteration
    p_vec(ib.pre) = pnew(ib.pre);
    
    % 6. Evaluate continuity equation
    divVel =  norm(Dx*unew + Dy*vnew);
    res(1,t) = divVel;
    
    % 7. Check for convergence --> start loop again or break
    if divVel < time.tol
        message = "Tolerance reached, stopping solver";
        disp(message);
        break;
    elseif divVel > 1e3
        message = "Solution diverging, stopping solver";
        disp(message);
        break;
    elseif t==time.iter
        message = "Maximum number of iterations reached";
        disp(message);
        break;
    else
        message = "No convergence reached in iteration number " + num2str(t) + ".";
        message2 = "Current residual is: " + num2str(divVel) + ".";
        message3 = "Continue with next itertion";
        disp(message);
        disp(message2);
        disp(message3);
        disp("---------------------------------------------");
    end
end

% Temporary plotting can be activated for debugging, will show flow
% fields, derivates and intermediate values

% subplot(3,3,1)
% surf(space.X,space.Y,reshape(u_vec,[space.dimY,space.dimX]))
% title('Velocity u')
% subplot(3,3,2)
% surf(space.X,space.Y,reshape(v_vec,[space.dimY,space.dimX]))
% title('Velocity v')
% subplot(3,3,3)
% surf(space.X,space.Y,reshape(p_vec,[space.dimY,space.dimX]))
% title('Pressure p')
% subplot(3,3,4)
% surf(space.X,space.Y,reshape(ustar,[space.dimY,space.dimX]))
% title('Intermediate Velocity u*')
% subplot(3,3,5)
% surf(space.X,space.Y,reshape(vstar,[space.dimY,space.dimX]))
% title('Intermediate Velocity v*')
% subplot(3,3,6)
% surf(space.X,space.Y,reshape(pc,[space.dimY,space.dimX]))
% title('Pressure correction pc')
% subplot(3,3,7)
% surf(space.X,space.Y,reshape(Dxustar,[space.dimY,space.dimX]))
% title('Dx*ustar with RC')
% subplot(3,3,8)
% surf(space.X,space.Y,reshape(Dyvstar,[space.dimY,space.dimX]))
% title('Dy*vstar with RC')
% subplot(3,3,9)
% surf(space.X,space.Y,reshape(Dxpc,[space.dimY,space.dimX]))
% title('Dx*pc');
% sgtitle('Variable Monitor')

% Reshape vectors for output to main script
u = reshape(u_vec,[space.dimY,space.dimX]);
v = reshape(v_vec,[space.dimY,space.dimX]);
p = reshape(p_vec,[space.dimY,space.dimX]);

end