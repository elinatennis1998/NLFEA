%Elina Geut 
%6/27/2019
% Compute Displacement and Traction Jumps for MRDG case

load('Node_U_V','Node_U_V');
load('GrainEps','GrainEps');
load('GrainDisp','GrainDisp');
load('GrainVol','GrainVol'); 

%Calculating Grain average based on Dr. Truster's notes
U_g(1,:) = GrainDisp(1,:) + GrainEps(1,:).*GrainXmid(1,:);
U_g(2,:) = GrainDisp(2,:) + GrainEps(2,:).*GrainXmid(2,:);

%Solid elements are arranged into Nodal array
Ug_new = zeros(numnp,2);
for i = 1:numelCG
        for j = 1:nnz(NodesOnElement(i,:))
            arrayCG = NodesOnElement(i,j);
            Ug_new(arrayCG,1) = U_g(1,i);
            Ug_new(arrayCG,2) = U_g(2,i);
        end
end

%DG elements are also arranged 
for i = 1:numel-numelCG
    condition1 = intersect(NodesOnElement(numelCG+i,1:4),NodesOnElement(1:numelCG,1:4));
    solid1 = find(any(sort(NodesOnElement(1:numelCG,1:nen_bulk),2) == (condition1'),2));
    for m = 1:3
        arrayDG = NodesOnElement(numelCG+i,2*nen_bulk+m);
        Ug_new(arrayDG,1) = U_g(1,solid1);
        Ug_new(arrayDG,2) = U_g(2,solid1);
    end
    
end

for i = 1:numel-numelCG
    condition2 = intersect(NodesOnElement(numelCG+i,5:8),NodesOnElement(1:numelCG,1:4));
    solid2 = find(any(sort(NodesOnElement(1:numelCG,1:nen_bulk),2) == (condition2'),2));
    for m = 1:3
        arrayDG = NodesOnElement(numelCG+i,2*nen_bulk+3+m);
        Ug_new(arrayDG,1) = U_g(1,solid2);
        Ug_new(arrayDG,2) = U_g(2,solid2);
    end
end
%Subtract U_g from the total nodal solution Node_U_V to get displacement
%jumps. 
Ug_new;
