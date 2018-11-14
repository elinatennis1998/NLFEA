function [NodeItem,numItem] = updateCGDGnodelist2(maxel,epnum,NodeCGDG,NodeMat,NodeItem_in,numItem_in)
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
    
    mats = find(NodeMat(node,:)>0);
    nodes = NodeMat(node,mats);
    numnodes = length(nodes);
    
    if numnodes == 1
        % Do only once
        for ncopy = 1:epnum(node)

            numItem = numItem + 1;
            nnode = NodeCGDG(ncopy,node);
            NodeItem(numItem,:) = row;
            NodeItem(numItem,1) = nnode;

        end
    
    else
        % Copy into next materials as well
        for i = 1:numnodes
            node2 = nodes(i); % get node ID in neighboring material
            for ncopy = 1:epnum(node2)

                numItem = numItem + 1;
                nnode = NodeCGDG(ncopy,node2);
                NodeItem(numItem,:) = row;
                NodeItem(numItem,1) = nnode;

            end
        end
        
    end
    
end

NodeItem = NodeItem(1:numItem,:);