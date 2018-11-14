function [SurfacesI,NodesOnElementt] = getixt(SurfacesI,numCL,NodesOnElement,nel)
%
% Tim Truster
% 05/01/12
% Function to generate NodesOnElementt array for 3D DG formulation, which is the table
% of integration triangles for the surfaces. The mesh is assumed to be
% "conforming" but discontinuous. Either bricks or tets of a single order
% are permitted.

% Assumes nel is same for all elements
if nel == 4 || nel == 10
    NodesOnElementt = zeros(numCL,3);
    SurI = zeros(numCL,2);
else
    NodesOnElementt = zeros(2*numCL,3);
    SurI = zeros(numCL,3);
end

for inter = 1:numCL
    
    elemL = SurfacesI(inter,5);
    edgeL = SurfacesI(inter,7);
    
    nelL = nel;
    
    %Extract patch nodal coordinates
    ElemFlagL = zeros(nelL, 1);
    for k = 1:nelL
        node = NodesOnElement(elemL,k);
        ElemFlagL(k) = node;
    end
    
    %Reorder element nodes in order to integrate on bottom side
    
    if nelL == 4
        if edgeL == 1	
            ilist = [2 4 3 1];
        elseif edgeL == 2
            ilist = [4 1 3 2];
        elseif edgeL == 3
            ilist = [1 4 2 3];
        elseif edgeL == 4
            ilist = [1 2 3 4];
        end
    end
    if nelL == 8
        if edgeL == 1	
            ilist = [5 1 4 8 6 2 3 7];
        elseif edgeL == 2
            ilist = [2 6 7 3 1 5 8 4];
        elseif edgeL == 3
            ilist = [5 6 2 1 8 7 3 4];
        elseif edgeL == 4
            ilist = [4 3 7 8 1 2 6 5];
        elseif edgeL == 5
            ilist = [1 2 3 4 5 6 7 8];
        else % edge == 6
            ilist = [8 7 6 5 4 3 2 1];
        end
    end
    if nelL == 10
        if edgeL == 1	
            ilist = [2 4 3 1 9 10 6 5 8 7];
        elseif edgeL == 2
            ilist = [4 1 3 2 8 7 10 9 5 6];
        elseif edgeL == 3
            ilist = [1 4 2 3 8 9 5 7 10 6];
        elseif edgeL == 4
            ilist = [1 2 3 4 5 6 7 8 9 10];
        end 
    end
    if nelL == 27
        if edgeL == 1	
            ilist = [5 1 4 8 6 2 3 7 17 12 20 16 18 10 19 14 13 9 11 15 23 24 22 21 25 26 27];
        elseif edgeL == 2
            ilist = [2 6 7 3 1 5 8 4 18 14 19 10 17 16 20 12 9 13 15 11 24 23 21 22 25 26 27];
        elseif edgeL == 3
            ilist = [5 6 2 1 8 7 3 4 13 18 9 17 15 19 11 20 16 14 10 12 25 26 23 24 22 21 27];
        elseif edgeL == 4
            ilist = [4 3 7 8 1 2 6 5 11 19 15 20 9 18 13 17 12 10 14 16 26 25 23 24 21 22 27];
        elseif edgeL == 5
            ilist = (1:27);
        else % edge == 6
            ilist = [8 7 6 5 4 3 2 1 15 14 13 16 11 10 9 12 20 19 18 17 22 21 23 24 26 25 27];
        end
    end
    ElemFlagL = ElemFlagL(ilist);
    
    if nel == 4 || nel == 10
        NodesOnElementt(inter,:) = ElemFlagL([1 3 2]);
        SurI(inter,:) = [1 inter];
    else
        NodesOnElementt(2*inter-1,:) = ElemFlagL([3 2 1]);
        NodesOnElementt(2*inter,:) = ElemFlagL([4 3 1]);
        SurI(inter,:) = [2 2*inter-1 2*inter];
    end
    
end

SurfacesI = [SurfacesI(:,1:8) SurI];