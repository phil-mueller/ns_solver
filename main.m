%% Preliminary stuff
clear
clc
close all 

%% Initialize variables and domain
disp("Initializing")
Init;

%% Build matrice operators
disp("Building the matrix derivative operators")
[Dx,Dy,Area,Length] = buildDxDy(space);

%% Build Laplacian Operator
disp("Building the matrix laplacian operators")
[Lu,Bu] = buildL(FV,space,boundary,value,'u');
[Lv,Bv] = buildL(FV,space,boundary,value,'v');
[Lp,Bp] = buildL(FV,space,boundary,value,'p');

 %% Solve the Navier Stokes equations for specified Reynolds numbers
% Specify desired Reynolds number --> use crescendo approach for large Re
ReNr=[25];
for i =1:size(ReNr,2)
    % Initialize
    clc
    ustart = u;
    vstart = v;
    pstart = p;
    % Adapt viscosity to fit Reynolds number
    Re = ReNr(i);
    FV.mu = (space.h*value.u_in*FV.rho)/Re;
    FV.nu = FV.mu/FV.rho;
    disp("Current Reynolds number is " + num2str(Re))
    % If final Reynolds number is reached, increase the number of
    % iterations for better results.
    if Re==ReNr(end)
        time.iter=3000;
        disp("Current Reynolds number is final Reynolds number")
        disp("Increasing the maximum number of iterations to " + num2str(time.iter))
    end
    % Call the solver
    [u,v,p,res] = solveNS(ustart,vstart,pstart,Dx,Dy,Lu,Lv,Lp,Bu,Bv,Bp,FV,time,space,value,ib,Area,Length);
end

%% Plot results
visualize(u,v,p,space,time,res,FV);