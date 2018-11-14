% Tim Truster
% 05/16/2015
%
% Subdivide a mesh of wedges into CONFORMING tetrahedra, linear elements
% only; important to keep track of the diagonals.
% Based on how to divide in tetrahedra.pdf

NodesOnElementT = zeros(3*numel,4);
RegionOnElementT = zeros(3*numel,1);

new_elem = 0;
for elem = 1:numel
    
    nodes = NodesOnElement(elem,1:6);
    % reorder nodes so that the lowest is in the lower-left corner
    minnode = min(nodes);
    if nodes(1) == minnode
        indirection = [1 2 3 4 5 6];
    elseif nodes(2) == minnode
        indirection = [2 3 1 5 6 4];
    elseif nodes(3) == minnode
        indirection = [3 1 2 6 4 5];
    elseif nodes(4) == minnode
        indirection = [4 6 5 1 3 2];
    elseif nodes(5) == minnode
        indirection = [5 4 6 2 1 3];
    elseif nodes(6) == minnode
        indirection = [6 5 4 3 2 1];
    end
    
    % Form the tetrahedra that we know
    new_elem = new_elem + 1;
    NodesOnElementT(new_elem,1:4) = nodes(indirection([1 5 6 4]));
    RegionOnElementT(new_elem) = RegionOnElement(elem);
    
    % Now decide on the other two
    Vi2Vi6 = min(nodes(indirection(2)),nodes(indirection(6)));
    Vi3Vi5 = min(nodes(indirection(3)),nodes(indirection(5)));
    
    if Vi2Vi6 < Vi3Vi5
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:4) = nodes(indirection([1 2 3 6]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:4) = nodes(indirection([1 2 6 5]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
    else
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:4) = nodes(indirection([1 2 3 5]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:4) = nodes(indirection([1 5 3 6]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
    end
    
end

NodesOnElement = NodesOnElementT;
RegionOnElement = RegionOnElementT;
numel = 3*numel;
nen = 4;
nen1 = 5;