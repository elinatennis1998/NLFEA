% Tim Truster
% 05/16/2015
%
% Subdivide a mesh of wedges into CONFORMING tetrahedra, quadratic elements
% only; important to keep track of the diagonals.
% Based on how to divide in tetrahedra.pdf

NodesOnElementT = zeros(3*numel,10);
RegionOnElementT = zeros(3*numel,1);

new_elem = 0;
for elem = 1:numel
    
    nodes = NodesOnElement(elem,1:18);
    % reorder nodes so that the lowest is in the lower-left corner
    minnode = min(nodes(1:6));
    if nodes(1) == minnode
        indirection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    elseif nodes(2) == minnode
        indirection = [2 3 1 5 6 4 8 9 7 11 12 10 14 15 13 17 18 16];
    elseif nodes(3) == minnode
        indirection = [3 1 2 6 4 5 9 7 8 12 10 11 15 13 14 18 16 17];
    elseif nodes(4) == minnode
        indirection = [4 6 5 1 3 2 12 11 10 9 8 7 13 15 14 18 17 16];
    elseif nodes(5) == minnode
        indirection = [5 4 6 2 1 3 10 12 11 7 9 8 14 13 15 16 18 17];
    elseif nodes(6) == minnode
        indirection = [6 5 4 3 2 1 11 10 12 8 7 9 15 14 13 17 16 18];
    end
    
    % Form the tetrahedra that we know
    new_elem = new_elem + 1;
    NodesOnElementT(new_elem,1:10) = nodes(indirection([1 5 6 4 16 11 18 13 10 12]));
    RegionOnElementT(new_elem) = RegionOnElement(elem);
    
    % Now decide on the other two
    Vi2Vi6 = min(nodes(indirection(2)),nodes(indirection(6)));
    Vi3Vi5 = min(nodes(indirection(3)),nodes(indirection(5)));
    
    if Vi2Vi6 < Vi3Vi5
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:10) = nodes(indirection([1 2 3 6 7 8 9 18 17 15]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:10) = nodes(indirection([1 2 6 5 7 17 18 16 14 11]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
    else
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:10) = nodes(indirection([1 2 3 5 7 8 9 16 14 17]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
        new_elem = new_elem + 1;
        NodesOnElementT(new_elem,1:10) = nodes(indirection([1 5 3 6 16 17 9 18 11 15]));
        RegionOnElementT(new_elem) = RegionOnElement(elem);
    end
    
end

NodesOnElement = NodesOnElementT;
RegionOnElement = RegionOnElementT;
numel = 3*numel;
nen = 10;
nen1 = 11;