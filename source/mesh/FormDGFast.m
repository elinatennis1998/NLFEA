function [ix,ixDG,numelDG,nummat,MatTypeTable,MateT,nen,nen1] = ...
         FormDGFast(SurfacesI,numCL,nen_bulk,mateprop,iel,ndm,maxmat,nonlin, ...
                ix,ixDG,numelDG,nummat,MatTypeTable,MateT)
%
% Tim Truster
% 08/07/2015
%
% Convert SurfaceI elements into DG elements for use in MATLAB or FEAP
% Combines the two elements from left and right sides into a single element
% with expanded ix array entry.
% Each combination of materials and nels for the elements on either side is
% treated as a new material set that is appended to the mesh.
% Code can be called multiple times to add more sets of interfaces; each
% interface input should be for a single type of interface element iel and
% a single set of material properties mateprop for the interface element.
%
% nen_bulk    = the max # of nodes per element on a bulk domain element
% maxmat      = maximum # of different DG element combinations to insert 
%               (e.g. (nelL=3,nelR=3); (nelL=4,nelR=3) => maxmat = 2)
%
% Faster version, keeps ixDG separate from the original ix, and no
% transpose is performed. Thus, ixDG = zeros(2*nen_bulk+1,0) needs to be
% initialized.

% Adjust size of ix, nen, MateT, nummat, ElemMatTable
% ix = ix';
nen1_old = size(ix,2);
RegionOnElement = ix(:,nen1_old);
nen = 2*nen_bulk;
nen1 = nen + 1;
ixDG = [ixDG zeros(nen1,numCL)];
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

for inter = 1:numCL
    
    elem1 = SurfacesI(inter,1);
    elem2 = SurfacesI(inter,2);
    fac1 = SurfacesI(inter,3);
    fac2 = SurfacesI(inter,4);
    
    nel1 = nnz(ix(elem1,1:nen_bulk));
    
    nel2 = nnz(ix(elem2,1:nen_bulk));
    
    %Extract patch nodal coordinates
    CouplerSet1 = ix(elem1,1:nel1);
    
    %Reorder element nodes in order to integrate on bottom side
    if ndm == 2
    if fac1 > 1
        if nel1 == 4
            if fac1 == 2
                ilist = [2 3 4 1];
            elseif fac1 == 3
                ilist = [3 4 1 2];
            else %edge == 4
                ilist = [4 1 2 3];
            end
        elseif nel1 == 3
            if fac1 == 2
                ilist = [2 3 1];
            elseif fac1 == 3
                ilist = [3 1 2];
            end
        elseif nel1 == 9
            if fac1 == 2
                ilist = [2 3 4 1 6 7 8 5 9];
            elseif fac1 == 3
                ilist = [3 4 1 2 7 8 5 6 9];
            else %edge == 4
                ilist = [4 1 2 3 8 5 6 7 9];
            end
        elseif nel1 == 6
            if fac1 == 2
                ilist = [2 3 1 5 6 4];
            elseif fac1 == 3
                ilist = [3 1 2 6 4 5];
            end
        end
        CouplerSet1 = CouplerSet1(ilist);
    end
    elseif ndm == 3
        if nel1 == 4
            if fac1 == 1	
                ilist = [2 4 3 1];
            elseif fac1 == 2
                ilist = [4 1 3 2];
            elseif fac1 == 3
                ilist = [1 4 2 3];
            elseif fac1 == 4
                ilist = [1 2 3 4];
            end
        end
        if nel1 == 8
            if fac1 == 1	
                ilist = [5 1 4 8 6 2 3 7];
            elseif fac1 == 2
                ilist = [2 6 7 3 1 5 8 4];
            elseif fac1 == 3
                ilist = [5 6 2 1 8 7 3 4];
            elseif fac1 == 4
                ilist = [4 3 7 8 1 2 6 5];
            elseif fac1 == 5
                ilist = [1 2 3 4 5 6 7 8];
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1];
            end
        end
        if nel1 == 10
            if fac1 == 1	
                ilist = [2 4 3 1 9 10 6 5 8 7];
            elseif fac1 == 2
                ilist = [4 1 3 2 8 7 10 9 5 6];
            elseif fac1 == 3
                ilist = [1 4 2 3 8 9 5 7 10 6];
            elseif fac1 == 4
                ilist = [1 2 3 4 5 6 7 8 9 10];
            end 
        end
        if nel1 == 27
            if fac1 == 1	
                ilist = [5 1 4 8 6 2 3 7 17 12 20 16 18 10 19 14 13 9 11 15 23 24 22 21 25 26 27];
            elseif fac1 == 2
                ilist = [2 6 7 3 1 5 8 4 18 14 19 10 17 16 20 12 9 13 15 11 24 23 21 22 25 26 27];
            elseif fac1 == 3
                ilist = [5 6 2 1 8 7 3 4 13 18 9 17 15 19 11 20 16 14 10 12 25 26 23 24 22 21 27];
            elseif fac1 == 4
                ilist = [4 3 7 8 1 2 6 5 11 19 15 20 9 18 13 17 12 10 14 16 26 25 23 24 21 22 27];
            elseif fac1 == 5
                ilist = (1:27);
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1 15 14 13 16 11 10 9 12 20 19 18 17 22 21 23 24 26 25 27];
            end
        end
        CouplerSet1 = CouplerSet1(ilist);
    end
    
    %Extract patch nodal coordinates
    CouplerSet2 = ix(elem2,1:nel2);
    
    %Reorder element nodes in order to integrate on bottom side
    if ndm == 2
    if fac2 > 1
        if nel2 == 4
            if fac2 == 2
                ilist = [2 3 4 1];
            elseif fac2 == 3
                ilist = [3 4 1 2];
            else %edge == 4
                ilist = [4 1 2 3];
            end
        elseif nel2 == 3
            if fac2 == 2
                ilist = [2 3 1];
            elseif fac2 == 3
                ilist = [3 1 2];
            end
        elseif nel2 == 9
            if fac2 == 2
                ilist = [2 3 4 1 6 7 8 5 9];
            elseif fac2 == 3
                ilist = [3 4 1 2 7 8 5 6 9];
            else %edge == 4
                ilist = [4 1 2 3 8 5 6 7 9];
            end
        elseif nel2 == 6
            if fac2 == 2
                ilist = [2 3 1 5 6 4];
            elseif fac2 == 3
                ilist = [3 1 2 6 4 5];
            end
        end
        CouplerSet2 = CouplerSet2(ilist);
    end
    elseif ndm == 3
        if nel2 == 4
            if fac2 == 1	
                ilist = [2 4 3 1];
            elseif fac2 == 2
                ilist = [4 1 3 2];
            elseif fac2 == 3
                ilist = [1 4 2 3];
            elseif fac2 == 4
                ilist = [1 2 3 4];
            end
        end
        if nel2 == 8
            if fac2 == 1	
                ilist = [5 1 4 8 6 2 3 7];
            elseif fac2 == 2
                ilist = [2 6 7 3 1 5 8 4];
            elseif fac2 == 3
                ilist = [5 6 2 1 8 7 3 4];
            elseif fac2 == 4
                ilist = [4 3 7 8 1 2 6 5];
            elseif fac2 == 5
                ilist = [1 2 3 4 5 6 7 8];
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1];
            end
        end
        if nel2 == 10
            if fac2 == 1	
                ilist = [2 4 3 1 9 10 6 5 8 7];
            elseif fac2 == 2
                ilist = [4 1 3 2 8 7 10 9 5 6];
            elseif fac2 == 3
                ilist = [1 4 2 3 8 9 5 7 10 6];
            elseif fac2 == 4
                ilist = [1 2 3 4 5 6 7 8 9 10];
            end 
        end
        if nel2 == 27
            if fac2 == 1	
                ilist = [5 1 4 8 6 2 3 7 17 12 20 16 18 10 19 14 13 9 11 15 23 24 22 21 25 26 27];
            elseif fac2 == 2
                ilist = [2 6 7 3 1 5 8 4 18 14 19 10 17 16 20 12 9 13 15 11 24 23 21 22 25 26 27];
            elseif fac2 == 3
                ilist = [5 6 2 1 8 7 3 4 13 18 9 17 15 19 11 20 16 14 10 12 25 26 23 24 22 21 27];
            elseif fac2 == 4
                ilist = [4 3 7 8 1 2 6 5 11 19 15 20 9 18 13 17 12 10 14 16 26 25 23 24 21 22 27];
            elseif fac2 == 5
                ilist = (1:27);
            else % edge == 6
                ilist = [8 7 6 5 4 3 2 1 15 14 13 16 11 10 9 12 20 19 18 17 22 21 23 24 26 25 27];
            end
        end
        CouplerSet2 = CouplerSet2(ilist);
    end
    
    
%%%%%%%%%%%%%%%%% IN FUTURE, add the reorientation of element 2 to match
%%%%%%%%%%%%%%%%% the DG templates
    
    
    
    % Add element to connectivity
    numelDG = numelDG + 1;
    ixDG(1:nel1+nel2,numelDG) = [CouplerSet1'; CouplerSet2'];
    
    % Add DG element pair to material list
    DGpair = [RegionOnElement(elem1) RegionOnElement(elem2) nel1 nel2];
    [tf, index] = ismember(DGpair, DGList, 'rows');
    if tf % copy other 
        ixDG(nen1,numelDG) = nummat_old + index;
    else % add to list
        nummat = nummat + 1;
        DGList = [DGList; DGpair];
        MateT(1:4+sizemateprop,nummat) = [DGpair'; mateprop'];
        MatTypeTable(1:3,nummat) = [nummat; iel; nonlin];
        ixDG(nen1,numelDG) = nummat;
    end
    
end


% Output only the actual number of added DG material IDs, and transpose the
% arrays back to the expected NL_FEA_Program format
MatTypeTable = MatTypeTable(:,1:nummat);
MateT = MateT(1:ndd,1:nummat)';