
% Average displacement and strain for linear elements 
load('GrainEps','GrainEps');
load('GrainVol','GrainVol');
load('Node_U_V','Node_U_V'); 

%Gerring the displacement average for each grain (note: not using shape
%functions, since element types are linear
for k = 1:numelCG
    var = NodesOnElement(k,1:nnz(NodesOnElement(k,:)));
    for m = 1:length(var)
        Displ_x(m) = Node_U_V(var(m),1);
        Displ_y(m) = Node_U_V(var(m),2);
    end
    GrainDispA(k,1) = sum(Displ_x)/nen_bulk;
    GrainDispA(k,2) = sum(Displ_y)/nen_bulk;
end
