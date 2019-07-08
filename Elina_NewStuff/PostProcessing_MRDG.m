%Elina Geut 
%6/27/2019
% Compute Displacement and Traction Jumps for MRDG case

load('Node_U_V','Node_U_V');
load('GrainEps','GrainEps');
load('GrainDisp','GrainDisp');
load('GrainVol','GrainVol'); 

nan = find(any(isnan(Node_U_V),2));
for k = 1:length(nan)
    Node_U_V(nan(k),1:2) = 0; %Replacing all NaN inputs with zeros
end
%Calculating Grain average based on Dr. Truster's notes

grainG = zeros(numgrain,numelemg);
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
        GrainEps(1,grainG(g,:)) = GrainEps(1,g);
        GrainEps(2,grainG(g,:)) = GrainEps(2,g);
        GrainDisp(1,grainG(g,:)) = GrainDisp(1,g);
        GrainDisp(2,grainG(g,:)) = GrainDisp(2,g);
        GrainXmid(1,grainG(g,:)) = GrainXmid(1,g);
        GrainXmid(2,grainG(g,:)) = GrainXmid(2,g);
    end

CoordinatesF = Coordinates';
% GrainXnode = zeros(ndm,(size(NodesOnElement(1:numelCG,1:nen_bulk),2))*(size(NodesOnElement(1:numelCG,1:nen_bulk),1)));
% GrainCoord = zeros(ndm,(size(NodesOnElement(1:numelCG,1:nen_bulk),2))*(size(NodesOnElement(1:numelCG,1:nen_bulk),1)));
% GrainDispl = zeros(ndm,(size(NodesOnElement(1:numelCG,1:nen_bulk),2))*(size(NodesOnElement(1:numelCG,1:nen_bulk),1)));
% GrainEpsl = zeros(ndm,(size(NodesOnElement(1:numelCG,1:nen_bulk),2))*(size(NodesOnElement(1:numelCG,1:nen_bulk),1)));

for i = 1:numelCG
    for j = 1:nen_bulk
        array = NodesOnElement(i,j);
        GrainXnode(1,array) = GrainXmid(1,i);
        GrainXnode(2,array) = GrainXmid(2,i);
        GrainCoord(1,array) = CoordinatesF(1,array);
        GrainCoord(2,array) = CoordinatesF(2,array);
        GrainDispl(1,array) = GrainDisp(1,i);
        GrainDispl(2,array) = GrainDisp(2,i);
        GrainEpsl(1,array) = GrainEps(1,i);
        GrainEpsl(2,array) = GrainEps(2,i);
    end
end

GrainXref = GrainCoord - GrainXnode; %Reference coordinate of the grain for GrainEps
Ug_new = zeros(numnp,2);

U_gg = GrainDispl + GrainEpsl.*GrainXref; %Calcullate average grain displacement in x
Ug_new = U_gg';

%Solid elements are arranged into Nodal array
for i = 1:length(U_gg)
    Ug_new(i,1) = U_gg(1,i);
    Ug_new(i,2) = U_gg(2,i);
end

%DG elements are also arranged 
for i = 1:numel-numelCG
    Ug_new(NodesOnElement(numelCG+i,nen_bulk*2+1:nen_bulk*2+3),1) = 0;
    Ug_new(NodesOnElement(numelCG+i,nen_bulk*2+1:nen_bulk*2+3),2) = 0;
    condition3 = sort(NodesOnElement(numelCG+i,1:nen_bulk),2);
     if NodesOnElement(numelCG+i,1:nen_bulk) == NodesOnElement(numelCG+i,nen_bulk+1:nen_bulk*2)
         solid3 = find(sum(sort(NodesOnElement(1:numelCG,1:nen_bulk),2) == [condition3],2) >= nen_bulk);
         arrayDG = NodesOnElement(numelCG+i,2*nen_bulk+3+3);
         Ug_new(arrayDG,1) = 0;
         Ug_new(arrayDG,2) = 0;
         array1 = NodesOnElement(numelCG+i,2*nen_bulk+3+1);
         arrayDG1 = find(NodeBC == NodesOnElement(numelCG+i,2*nen_bulk+3+1));
         Ug_new(array1,1) = NodeBC(arrayDG1(1),3);
         Ug_new(array1,2) = NodeBC(arrayDG1(2),3);
         array2 = NodesOnElement(numelCG+i,2*nen_bulk+3+2);
         arrayDG2 = find(NodeBC == NodesOnElement(numelCG+i,2*nen_bulk+3+2));
         Ug_new(array2,1) = NodeBC(arrayDG2(1),3);
         Ug_new(array2,2) = NodeBC(arrayDG2(2),3);
         array3 = NodesOnElement(numelCG+i,2*nen_bulk+3+3);
         arrayDG3 = find(NodeBC == NodesOnElement(numelCG+i,2*nen_bulk+3+3));
         Ug_new(array3,1) = NodeBC(arrayDG3(1),3);
         Ug_new(array3,2) = NodeBC(arrayDG3(2),3);
    end
end

%Subtract U_g from the total nodal solution Node_U_V to get displacement
%jumps. 
Ug_new;
DispFine = Node_U_V - Ug_new;
%Save the results to impose onto CS
save('DispFine','DispFine');