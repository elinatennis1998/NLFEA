%Elina Geut
% Some essential plotting functions

plotElemCont2(Coordinates,StreListE(1,1:numelCG),NodesOnElement(1:numelCG,1:4),1,(1:size(NodesOnElement,1)),[1 0 0])
plotElemCont2(Coordinates,StreListCF(1,1:numelCG),NodesOnElement(1:numelCG,1:4),1,(1:size(NodesOnElement,1)),[1 0 0])
plotNodeCont2(Coordinates + Node_U_V*2, Node_U_V(:,1), NodesOnElement, 1, (1:numelCG))
plotMesh2(Coordinates,NodesOnElement,1,(1:size(NodesOnElement,1)),[1 1 1 1 0])