function [ix,numel,nummat,MatTypeTable,MateT,nen,nen1] = ...
         FormDGP1P0(SurfacesI,numCL,nen_bulk,mateprop,iel,ndm,maxmat,nonlin, ...
                ix,numel,nummat,MatTypeTable,MateT)
%
% Tim Truster
% 04/17/2013
%
% Convert SurfaceI elements into DG elements for use in MATLAB or FEAP
% Combines the two elements from left and right sides into a single element
% with expanded ix array entry.
% Each combination of materials and nels for the elements on either side is
% treated as a new material set that is appended to the mesh.
% Code can be called multiple times to add more sets of interfaces; each
% interface input should be for a single type of interface element iel and
% a single set of material properties mateprop for the interface element.

% nen_bulk    = the max # of nodes per element on a bulk domain element
% maxmat      = maximum # of different DG element combinations to insert 
%               (e.g. (nelL=3,nelR=3); (nelL=4,nelR=3) => maxmat = 2)

% Adjust size of ix, nen, MateT, nummat, ElemMatTable
ix = ix';
nen1_old = size(ix,1);
nen = 2*nen_bulk;
nen1 = nen + 1;
ix = [[ix(1:nen1_old-1,1:numel); zeros(nen1-nen1_old,numel); ix(nen1_old,1:numel)] zeros(nen1,numCL)];
MateT = MateT';
ndd_old = size(MateT,1);
sizemateprop = size(mateprop,2);
ndd_new = 4 + sizemateprop;
if ndd_new > ndd_old
    ndd = ndd_new;
    MateT = [[MateT; zeros(ndd_new-ndd_old,nummat)] zeros(ndd,maxmat)];
else
    ndd = ndd_old;
    MateT = [MateT zeros(ndd,maxmat)];
end
sizeMTT = size(MatTypeTable,1);
MatTypeTable = [MatTypeTable zeros(sizeMTT,maxmat)];
nummat_old = nummat;

DGList = zeros(0,4);

if ndm == 2 || ndm == 3

for inter = 1:numCL
    
    elemL = SurfacesI(inter,1);
    elemR = SurfacesI(inter,2);
    edgeL = SurfacesI(inter,3);
    edgeR = SurfacesI(inter,4);
    
    nelL = nnz(ix(1:nen,elemL));
%     nelL = getnel(ix(:,elemL),nen_bulk,ndm);
    
    nelR = nnz(ix(1:nen,elemR));
%     nelR = getnel(ix(:,elemR),nen_bulk,ndm);
    
    %Extract patch nodal coordinates
    ElemFlagL = ix(1:nelL,elemL);
    
    %Reorder element nodes in order to integrate on bottom side
    if ndm == 2
    if edgeL > 1
        if nelL == 4
            if edgeL == 2
                ilist = [2 3 1 4];
            elseif edgeL == 3
                ilist = [3 1 2 4];
            end
        elseif nelL == 3
            if edgeL == 2
                ilist = [2 3 1];
            elseif edgeL == 3
                ilist = [3 1 2];
            end
        elseif nelL == 9
            if edgeL == 2
                ilist = [2 3 4 1 6 7 8 5 9];
            elseif edgeL == 3
                ilist = [3 4 1 2 7 8 5 6 9];
            else %edge == 4
                ilist = [4 1 2 3 8 5 6 7 9];
            end
        elseif nelL == 6
            if edgeL == 2
                ilist = [2 3 1 5 6 4];
            elseif edgeL == 3
                ilist = [3 1 2 6 4 5];
            end
        end
        ElemFlagL = ElemFlagL(ilist);
    end
    elseif ndm == 3
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
    end
    
    %Extract patch nodal coordinates
    ElemFlagR = ix(1:nelR,elemR);
    
    %Reorder element nodes in order to integrate on bottom side
    if ndm == 2
    if edgeR > 1
        if nelR == 4
            if edgeR == 2
                ilist = [2 3 1 4];
            elseif edgeR == 3
                ilist = [3 1 2 4];
            end
        elseif nelR == 3
            if edgeR == 2
                ilist = [2 3 1];
            elseif edgeR == 3
                ilist = [3 1 2];
            end
        elseif nelR == 9
            if edgeR == 2
                ilist = [2 3 4 1 6 7 8 5 9];
            elseif edgeR == 3
                ilist = [3 4 1 2 7 8 5 6 9];
            else %edge == 4
                ilist = [4 1 2 3 8 5 6 7 9];
            end
        elseif nelR == 6
            if edgeR == 2
                ilist = [2 3 1 5 6 4];
            elseif edgeR == 3
                ilist = [3 1 2 6 4 5];
            end
        end
        ElemFlagR = ElemFlagR(ilist);
    end
    elseif ndm == 3
        if nelR == 4
            if edgeR == 1	
                ilist = [2 4 3 1];
            elseif edgeR == 2
                ilist = [4 1 3 2];
            elseif edgeR == 3
                ilist = [1 4 2 3];
            elseif edgeR == 4
                ilist = [1 2 3 4];
            end
        end
        if nelR == 8
            if edgeR == 1	
                ilist = [5 1 4 8 6 2 3 7];
            elseif edgeR == 2
                ilist = [2 6 7 3 1 5 8 4];
            elseif edgeR == 3
                ilist = [5 6 2 1 8 7 3 4];
            elseif edgeR == 4
                ilist = [4 3 7 8 1 2 6 5];
            elseif edgeR == 5
                ilist = [1 2 3 4 5 6 7 8];
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1];
            end
        end
        if nelR == 10
            if edgeR == 1	
                ilist = [2 4 3 1 9 10 6 5 8 7];
            elseif edgeR == 2
                ilist = [4 1 3 2 8 7 10 9 5 6];
            elseif edgeR == 3
                ilist = [1 4 2 3 8 9 5 7 10 6];
            elseif edgeR == 4
                ilist = [1 2 3 4 5 6 7 8 9 10];
            end 
        end
        if nelR == 27
            if edgeR == 1	
                ilist = [5 1 4 8 6 2 3 7 17 12 20 16 18 10 19 14 13 9 11 15 23 24 22 21 25 26 27];
            elseif edgeR == 2
                ilist = [2 6 7 3 1 5 8 4 18 14 19 10 17 16 20 12 9 13 15 11 24 23 21 22 25 26 27];
            elseif edgeR == 3
                ilist = [5 6 2 1 8 7 3 4 13 18 9 17 15 19 11 20 16 14 10 12 25 26 23 24 22 21 27];
            elseif edgeR == 4
                ilist = [4 3 7 8 1 2 6 5 11 19 15 20 9 18 13 17 12 10 14 16 26 25 23 24 21 22 27];
            elseif edgeR == 5
                ilist = (1:27);
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1 15 14 13 16 11 10 9 12 20 19 18 17 22 21 23 24 26 25 27];
            end
        end
        ElemFlagR = ElemFlagR(ilist);
    end
    
    % Add element to connectivity
    numel = numel + 1;
    ix(1:nelL+nelR,numel) = [ElemFlagL; ElemFlagR];
    
    % Add DG element pair to material list
    DGpair = [ix(nen1,elemL) ix(nen1,elemR) nelL nelR];
    [tf, index] = ismember(DGpair, DGList, 'rows');
    if tf % copy other 
        ix(nen1,numel) = nummat_old + index;
    else % add to list
        nummat = nummat + 1;
        DGList = [DGList; DGpair];
        MateT(1:4+sizemateprop,nummat) = [DGpair'; mateprop'];
        MatTypeTable(1:3,nummat) = [nummat; iel; nonlin];
        ix(nen1,numel) = nummat;
    end
    
end

else
    
    disp('ndm ~= 2 or 3 is not allowed')
    return
    
end

% Output only the actual number of added DG material IDs, and transpose the
% arrays back to the expected NL_FEA_Program format
ix = ix';
MatTypeTable = MatTypeTable(:,1:nummat);
MateT = MateT(1:ndd,1:nummat)';