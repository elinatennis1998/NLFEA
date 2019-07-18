function [RegionOnElement,MatTypeTable,MateT,nummat,grainG,BoundGrains,CornerGrain] = GeomProp(numgrain,numelemg,tfact,numgs,numgh,bCrys,m2,m1,nel,nu)
%Elina Geut 6/27/2019
%   Function for assigning RVE properties 
grainG = zeros(numgrain,numelemg*tfact);
grain = 0;
% for each grain, which elements belong to it
for j = 1:numgs
    for i = 1:numgh
        el = 0;
        grain = grain + 1;
        for m = 1:bCrys
            for l = 1:bCrys*tfact
                if nel == 3
                    elem = (j-1)*nu*bCrys*(j*tfact-2*(j-1))+(i-1)*bCrys*tfact; % bottom-corner element of grain
                    elem = elem + (m-1)*nu*tfact;
                    elem = elem + (l-1) + 1;
                elseif nel == 4
                    elem = (j-1)*nu*bCrys+(i-1)*bCrys;
                    elem = elem + (m-1)*nu*tfact;
                    elem = elem + (l-1) + 1;
                end
                el = el + 1;
                grainG(grain,el) = elem;
            end
        end
    end
end
% % inverse map: the grain that an element belongs to
    for g = 1:numgrain
        RegionOnElement(grainG(g,:)) = g;
    end
    RegionOnElement = RegionOnElement';
%% Set up phase pattern and material properties
MatTypeTable = [1:numgrain; ones(1,numgrain)];
mats = [m1; m2];
alterphase = 2*rem(1:numgrain,2) - 1;
alterphase(alterphase==-1) = 2;
MateT = mats(alterphase,:);%% Output quantity flags
nummat = numgrain;
BoundGrains = zeros(numgh*numgs-(numgh-2)*(numgs-2),1);
corner = [1 numgh];
for m = 1:numgs
    for k = 1:numgh
        if m == 1
            BoundGrains(k) = k;
        elseif m == numgs
            BoundGrains(m*numgh-(numgh-k)) = numgs*(numgh-1)+k;
        else
            for s = 1:2
                BoundGrains(numgh*(m-1)+s) = corner(s)+numgh*(m-1);
            end
        end
    end
end

CornerGrain(1) = 1;
CornerGrain(2) = numgh;
CornerGrain(3) = numgh*(numgs-1)+1;
CornerGrain(4) = numgh*numgs;
end

