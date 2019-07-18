function [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel] = Meso_Locking_Truster(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,grainG,nen_bulk,numelemg,GrainA,numel,numnpMicro,numnpMeso,MateT,numgrain,nummat,MatTypeTable,RegionOnElement)
%Elina Geut
%Created 3/4/2019
%Last Modified 3/4/2019
%Reassigning BC with locked grains
 
for i = 1:num_locked_g
     
    micro = unique(reshape(NodesOnElement(grainG(locked_g(i),:),1:nen_bulk),numelemg*nen_bulk,1)); %Corrected
    l(i) = length(micro);
    r_g = numgrain-locked_g(i);
    NodeBC = [NodeBC
%         [(numnpMicro+1:numnpMicro+(locked_g(i)-1)*meso_nen)' ones((locked_g(i)-1)*meso_nen,1) zeros((locked_g(i)-1)*meso_nen,1)
%         (numnpMicro+1:numnpMicro+(locked_g(i)-1)*meso_nen)' 2*ones((locked_g(i)-1)*meso_nen,1) zeros((locked_g(i)-1)*meso_nen,1)]
%         [(numnpMicro+locked_g(i)*meso_nen+1:numnpMicro+numnpMeso)'   ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
%         (numnpMicro+locked_g(i)*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]
        [micro ones(l(i),1) zeros(l(i),1)
        micro 2*ones(l(i),1) zeros(l(i),1)]];
    % Correction for grain 5 of mesoscale
    index = find(MateT(:,1) == locked_g(i));
    MateT(index(1),5:8) = [1 0 0 1];
    MateT(index(2),5:8) = [1 0 0 1];
    index2 = find(MateT(:,2) == locked_g(i));
    MateT(index2(1),5:8) = [0 1 1 0];
    MateT(index2(2),5:8) = [0 1 1 0];
    l_MT = length(MateT)+1;
    MateT(l_MT,1:8) = [MateT(locked_g(i),1:3) GrainA 0 0 0 0];
    MatTypeTable(1:3,l_MT) = [l_MT 11 0]';
    numel = numel + 1;
    %Recised up until this point
    NodesOnElement(numel,1:12) = [numnpMicro+(locked_g(i)-1)*meso_nen+1 numnpMicro+(locked_g(i)-1)*meso_nen+2 numnpMicro+(locked_g(i)-1)*meso_nen+3 zeros(1,12-3)];
    RegionOnElement(numel) = l_MT;
    nummat = nummat + 1;
end
 
    r_g = numgrain-locked_g(2);
    NodeBC = [NodeBC
        [(numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)
        (numnpMicro+1:numnpMicro+(locked_g(1)-1)*meso_nen)' 2*ones((locked_g(1)-1)*meso_nen,1) zeros((locked_g(1)-1)*meso_nen,1)]
        [(numnpMicro+locked_g(1)*meso_nen+1:numnpMicro+(locked_g(2)-1)*meso_nen)'   ones(1*meso_nen,1) zeros(1*meso_nen,1)
        (numnpMicro+locked_g(1)*meso_nen+1:numnpMicro+(locked_g(2)-1)*meso_nen)' 2*ones(1*meso_nen,1) zeros(1*meso_nen,1)]
        [(numnpMicro+locked_g(2)*meso_nen+1:numnpMicro+numnpMeso)'   ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)
        (numnpMicro+locked_g(2)*meso_nen+1:numnpMicro+numnpMeso)' 2*ones(r_g*meso_nen,1) zeros(r_g*meso_nen,1)]];
    numBC = length(NodeBC);
 
end