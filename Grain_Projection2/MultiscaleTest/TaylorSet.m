% Tim Truster
% 06/30/2017
%
% Set up fluxes for Taylor model
isw = 54;
ndf = ndm;
psflags

% CoordinatesI

for interI = 1:MaterTypeNum(3)-MaterTypeNum(2)
    
    trjump = BoundTay(1:18,interI);
    if ndm == 2
        trjump = reshape(trjump(1:6),2,3);
    else
        error('nmd=3')
    end
    fluxT = factorT*trjump*eRVE';
    CoordinatesI(2*(ndm+1)*(interI-1)+1,:) = fluxT';
        
end