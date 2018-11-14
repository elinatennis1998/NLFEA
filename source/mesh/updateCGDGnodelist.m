function [NodeItem,numItem] = updateCGDGnodelist(maxel,epnum,NodeCGDG,NodeItem_in,numItem_in)
%
% Tim Truster
% 01/01/2014
%
% Convert a list of nodal entries (e.g. boundary conditions) on a CG mesh
% into a DG mesh created by CGtoDGmesh

NodeItem = zeros(maxel*numItem_in,size(NodeItem_in,2));
numItem = 0;

for entry = 1:numItem_in
    
    node = NodeItem_in(entry,1);
    row = NodeItem_in(entry,:);
    for ncopy = 1:epnum(node)
        
        numItem = numItem + 1;
        nnode = NodeCGDG(ncopy,node);
        NodeItem(numItem,:) = row;
        NodeItem(numItem,1) = nnode;
        
    end
    
end

NodeItem = NodeItem(1:numItem,:);