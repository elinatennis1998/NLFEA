% Tim Truster
% 10/13/2013
%
% Continuation for Elasticity RFB calculation - for bubble in y direction
% Called by RFBEfuncE

newloads = 1;
numSL_new = numSL;
SurfacesL_new = zeros(numSL_new,7);
SurfacesL_new(1,:) = [4 1 1 1 0 1 0];
for i = 2:numSL_new-1;
    SurfacesL_new(i,:) = [3+i 2+i 2*i-1 1 -1 1 0];
end
SurfacesL_new(numSL,:) = [2 2+minc 2*minc-1 1 -1 1 0];