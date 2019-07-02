%Elina Geut 
%6/27/2019
% Compute Displacement and Traction Jumps for MRDG case

load('Node_U_V','Node_U_V');
load('GrainEps','GrainEps');
load('GrainDisp','GrainDisp');
load('GrainVol','GrainVol'); 

nan = find(any(isnan(Node_U_V),2));
for k = 1:length(nan)
Node_U_V(nan(k),1:2) = 0;
end
%Calculating Grain average based on Dr. Truster's notes
U_g(1,:) = GrainDisp(1,:) + GrainEps(1,:).*GrainXmid(1,:); %Calcullate average grain displacement in x
U_g(2,:) = GrainDisp(2,:) + GrainEps(2,:).*GrainXmid(2,:); %Calcullate average grain displacement in y

%Solid elements are arranged into Nodal array
Ug_new = zeros(numnp,2);
U_gg = zeros(ndf,numgrain*2);

for i = 1:numelCG/tfact
    U_gg(1,i*tfact-1:i*tfact) = U_g(1,i);
    U_gg(2,i*tfact-1:i*tfact) = U_g(2,i);
end
for i = 1:numelCG
    for j = 1:nnz(NodesOnElement(i,:))
        arrayCG = NodesOnElement(i,j);
        Ug_new(arrayCG,1) = U_gg(1,i);
        Ug_new(arrayCG,2) = U_gg(2,i);
    end
end

%DG elements are also arranged 
for i = 1:numel-numelCG
    condition1 = intersect(NodesOnElement(numelCG+i,1:nen_bulk),NodesOnElement(1:numelCG,1:nen_bulk));
    solid1 = find(sum(sort(NodesOnElement(1:numelCG,1:nen_bulk),2) == [condition1]',2)==nen_bulk);
    for m = 1:3
        arrayDG = NodesOnElement(numelCG+i,2*nen_bulk+m);
        Ug_new(arrayDG,1) = U_gg(1,solid1);
        Ug_new(arrayDG,2) = U_gg(2,solid1);
    end
    
end

for i = 1:numel-numelCG
    condition2 = intersect(NodesOnElement(numelCG+i,nen_bulk+1:nen_bulk+3),NodesOnElement(1:numelCG,1:nen_bulk));
    solid2 = find(sum(sort(NodesOnElement(1:numelCG,1:nen_bulk),2) == [condition2]',2)==nen_bulk);
    for m = 1:3
        arrayDG = NodesOnElement(numelCG+i,2*nen_bulk+3+m);
        Ug_new(arrayDG,1) = U_gg(1,solid2);
        Ug_new(arrayDG,2) = U_gg(2,solid2);
    end
end
%Subtract U_g from the total nodal solution Node_U_V to get displacement
%jumps. 
Ug_new;
DispFine = Node_U_V - Ug_new;
%Save the results to impose onto CS
save('DispFine','DispFine');