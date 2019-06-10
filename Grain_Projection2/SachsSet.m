% Tim Truster
% 06/30/2017
%
% Set up fluxes for Taylor model
isw = 54;
ndf = ndm;
psflags

% CoordinatesI

for interI = 1:MaterTypeNum(3)-MaterTypeNum(2)
    
    matLR = RegionsOnInterface(interI,2:3);
    trjumpi = BoundSachs(1:18,interI);
    if ndm == 2
        iDmatL = reshape(trjumpi(1:9),3,3);
        iDmatR = reshape(trjumpi(10:18),3,3);
    else
        error('nmd=3')
    end
    xint = GrainXmid(1:2,matLR(1)) - GrainXmid(1:2,matLR(2));
    xvect = [xint(1) 0 xint(2)/2
             0 xint(2) xint(1)/2];
%     jumpS = (xvect + trjumpi*DmatSachs)*eRVE';
%     CoordinatesI((ndm+2)*(interI-1)+2,:) = jumpS';
    jumpSu = factorS*xvect*eRVE';
    jumpSe = factorS*DmatSachs*eRVE';
    CoordinatesI(2*(ndm+1)*(interI-1)+4,:) = jumpSu';
    CoordinatesI(2*(ndm+1)*(interI-1)+5,:) = jumpSe(1:2)';
    CoordinatesI(2*(ndm+1)*(interI-1)+6,:) = [jumpSe(3) 0];
        
end