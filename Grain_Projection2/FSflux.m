% Tim Truster
% 07/07/2017
%
% Compute fine scale fluxes
GrainInteg
InterInteg

for interI = MaterTypeNum(2)-nummatCG:MaterTypeNum(3)-nummatCG
    
    CoordinatesI(2*(ndm+1)*(interI-1)+1,:) = BoundSig(5:6,interI)' - BoundSig(9:10,interI)'; %fluxjc
    CoordinatesI(2*(ndm+1)*(interI-1)+2,:) = BoundSig(23:24,interI)'; %fluxjx
    CoordinatesI(2*(ndm+1)*(interI-1)+3,:) = BoundSig(25:26,interI)'; %fluxjy
    CoordinatesI(2*(ndm+1)*(interI-1)+4,:) = BoundSig(3:4,interI)'; %jumpiu
    CoordinatesI(2*(ndm+1)*(interI-1)+5,:) = BoundSig(19:20,interI)'; %jumpie part 1
    CoordinatesI(2*(ndm+1)*(interI-1)+6,:) = BoundSig(21:22,interI)'; %jumpie part 2
    
end
save('CoordinatesI.mat','CoordinatesI');