%Elina Geut 
%Script for imposing the FS stresses onto the CS
%Created 5/10/2019
%Last modified 5/10/2019

count = tfact*bCrys^2*num_locked_g;
StreListCF = StreListE;

for j = 1:num_locked_g
    grain_e = find(RegionOnElement (:,1) == locked_g(j));
    for i = 1:count
        StreListCF(1,grain_e) = StreListE(1,ccell(j));
%         StreListCF(1,(grain_e1:grain_e2)) = StreListE(1,ccell(j));
    end
end
 %plotElemCont2(Coordinates,StreListCF(1,1:numelCG),NodesOnElement(1:numelCG,1:4),1,(1:size(NodesOnElement,1)),[1 0 0])