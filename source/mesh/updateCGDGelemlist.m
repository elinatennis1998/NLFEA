function [SurfaceItem,numItem] = updateCGDGelemlist(ndm,nen,ixCG,ixDG,SurfaceItem_in,numItem_in)
%
% Tim Truster
% 01/01/2014
%
% Convert a list of surface entries (e.g. tractions) on a CG mesh
% into a DG mesh created by CGtoDGmesh

SurfaceItem = SurfaceItem_in;
numItem = numItem_in;

if ndm == 2 % Update node IDs

for entry = 1:numItem_in
    
    elem = SurfaceItem_in(entry,3);
    
    for col = 1:2
        
        nodeCG = SurfaceItem_in(entry,col);
%         ind = find(ixCG(elem,1:nen)==nodeCG);
        nodeDG = ixDG(elem,ixCG(elem,1:nen)==nodeCG);%ind);
        SurfaceItem(entry,col) = nodeDG;
        
    end
    
end

end
