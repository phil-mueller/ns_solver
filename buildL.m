function [L,B,ib] = buildL(FV,space,boundary,value,str)
% This function builds the L matrix operator
% To build L, the finite volume laplacian formulation that is already
% available form sessions 4-6 of the course is used

% If the script is being called for the first time generate stencil files
if isfile('build_inner.m')
else
    generate_stencil_innernode;
    generate_stencil_north;
    generate_stencil_east;
    generate_stencil_south;
    generate_stencil_west;
    generate_stencil_northEast;
    generate_stencil_southEast;
    generate_stencil_southWest;
    generate_stencil_northWest;
end

% Index function and initialization
index = @(ii, jj) ii + (jj-1) * size(space.X, 1);
L = zeros(size(space.X,1)*size(space.X,2));
B = zeros(size(space.X));

% Loop over all nodes
for i = 1:size(space.X,1)
    for j = 1:size(space.X,2)
        [L(index(i,j),:),b] = stamp(i,j,space,FV,boundary,value,str);
        B(i,j) = b;
    end
end

% Make matrices sparse so save space
L = sparse(L);
B = sparse(B);

end

