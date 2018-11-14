function [xs, ecn, cis, ixs, numnps, cell] = AttachPatches(xs, ecn, cis, numP1, xc, el2mesh, ixs, ixc, ix, PFTable, numA, minc, numnps, cell, nel, cel, celn)
% [xs, ecnout, cisout, PatchConnout, Uout, Vout,
% Wout, ShelfTableout, ShelfIndexout, numPout] = AttachPatches(xs1, 
% ecn1, cis1, PatchConn1, U1, V1, W1, ShelfTable1, ShelfIndex1, 
% numP1, xc, ecn2, cis2, PatchConn2, U2, V2, W2, 
% ShelfTable2, ShelfIndex2, numP2, PFTable, numA)
%
% Tim Truster
% CEE Graduate Student
% UIUC
% 1/21/2009
%
%
%   This subroutine attaches a two models together (1 and 2). This means
%   that the Knot Vectors of the two models are combined, the Node Tables
%   are combined while also eliminating shared nodes at attachment faces,
%   and patch connectivities are updated to reflect new node numbers.
%   Compatibilities of the specified faces are checked before models are
%   combined.
%   (1) - Model to be retained; model to be appended; defined as Global
%   (2) - Patch or set of patches to be appended to Model(1); defined as Sub
%

%     --------------VARIABLE DECLARATIONS--------------------------------

%   Coordinates, PatchTable, PatchIndex, PatchConn, U, V, W, numP -> see NURBSGUI
%   PFTable = Table of patches and faces in Global and Sub Models of
%             interfaces for attachement; organized as follows:
%             [Patch1 Face1 Patch2 Face2]
%   numA = length(PFTable) = total number of attachments to make

%     -------------------------------------------------------------------

%   Sort table of attachments into ascending order of patches in Global Model
PFTable = PFTable';
[PFTable(1:numA,1),sind] = sort(PFTable(1:numA,1));
for i = 2:4
    PFTable(1:numA,i) = PFTable(sind,i);
end
PFTable = PFTable';

if nel == 27
    minc = 2*minc;
end
ncel3 = (minc + 1)^3;
ncel2 = (minc + 1)^2;
% numP2 = 1;
% ecn2 = zeros(ncel3,1);
% for i = 1:ncel3
%     ecn2(i) = i;
% end
% cis2 = [1 nel
%         ncel3+1 nel];
TotalNodeList = zeros(2,ncel2*6+nel); %Aggregate list of all duplicated nodes
ecn2 = zeros(celn,1);
% NodeList = zeros(ncel2,2); %List of nodes that will be duplicates

for i = 1:nel
    TotalNodeList(1,i) = el2mesh(i);
    TotalNodeList(2,i) = ix(i);
end
lenNL = nel;

%Make numA attachments if possible
for att = 1:numA
	
    %Load patch and face data
    Patch1 = PFTable(1,att);
    Face1 = PFTable(2,att);
%     Patch2 = PFTable(3,att);
    Face2 = PFTable(4,att);
    
    %   Load parameter data for matching faces
% 	nu1 = minc + 1;
% 	nv1 = minc + 1;
% 	nw1 = minc + 1;
% 	nu2 = minc + 1;
% 	nv2 = minc + 1;
% 	nw2 = minc + 1;
    %   Identify parameter directions of patch in Global model on face of
    %   attachment
    if Face1 < 3 %Face = 1 or 2; V & W on face
%         n(1,3) = nu1;
%         n(1,1) = nw1;
%         n(1,2) = nv1;
        if Face1 == 1 %negative face
            ind13 = 1;
        else %positive face
            ind13 = minc + 1;
        end
    elseif Face1 > 4 %Face = 5 or 6; U & V on face
%         n(1,3) = nw1;
%         n(1,1) = nv1;
%         n(1,2) = nu1;
        if Face1 == 5 %negative face
            ind13 = 1;
        else %positive face
            ind13 = minc + 1;
        end
    else %Face = 3 or 4; U & W on face
%         n(1,3) = nv1;
%         n(1,1) = nw1;
%         n(1,2) = nu1;
        if Face1 == 3 %negative face
            ind13 = 1;
        else %positive face
            ind13 = minc + 1;
        end
    end
    %   Identify parameter directions of patch in Sub model on face of
    %   attachment
    if Face2 < 3 %Face = 1 or 2; V & W on face
%         n(2,3) = nu2;
%         n(2,1) = nw2;
%         n(2,2) = nv2;
        if Face2 == 1 %negative face
            ind23 = 1;
        else %positive face
            ind23 = minc + 1;
        end
    elseif Face2 > 4 %Face = 5 or 6; U & V on face
%         n(2,3) = nw2;
%         n(2,1) = nv2;
%         n(2,2) = nu2;
        if Face2 == 5 %negative face
            ind23 = 1;
        else %positive face
            ind23 = minc + 1;
        end
    else %Face = 3 or 4; U & W on face
%         n(2,3) = nv2;
%         n(2,1) = nw2;
%         n(2,2) = nu2;
        if Face2 == 3 %negative face
            ind23 = 1;
        else %positive face
            ind23 = minc + 1;
        end
    end

    %Global Model - Knot Vectors, loop counter increments
    if Face1 < 3 %Face = 1 or 2; V & W on face
        inc11 = ncel2;
        inc12 = (minc + 1);
        inc13 = 1;
    elseif Face1 > 4 %Face = 5 or 6; U & V on face
        inc11 = (minc + 1);
        inc12 = 1;
        inc13 = ncel2;
    else %Face = 3 or 4; U & W on face
        inc11 = ncel2;
        inc12 = 1;
        inc13 = (minc + 1);
    end
    
    %Global Model - Nodes
    ind = cis(1,Patch1);
    Corners1 = zeros(4,3);
    Corners2 = zeros(4,3);
%     for k = 1:(minc + 1)
%         ind1 = (k-1)*inc11 + ind + (ind13 - 1)*inc13;
%         for j = 1:(minc + 1)
%             node = ecn(ind1);
%             for l = 1:3
%                 FacePoints1(k,j,l) = xs(l,node);
%             end
%             ind1 = ind1 + inc12;
%         end
%     end
    ind1 = ind + (ind13 - 1)*inc13;
    node = ecn(ind1);
    for l = 1:3
        Corners1(1,l) = xs(l,node);
    end
    ind1 = ind1 + inc12*minc;
    node = ecn(ind1);
    for l = 1:3
        Corners1(2,l) = xs(l,node);
    end
    ind1 = minc*inc11 + ind + (ind13 - 1)*inc13;
    node = ecn(ind1);
    for l = 1:3
        Corners1(4,l) = xs(l,node);
    end
    ind1 = ind1 + inc12*minc;
    node = ecn(ind1);
    for l = 1:3
        Corners1(3,l) = xs(l,node);
    end

    %Sub Model - Knot Vectors, loop counter increments
    if Face2 < 3 %Face = 1 or 2; V & W on face
        inc21 = ncel2;
        inc22 = (minc + 1);
        inc23 = 1;
    elseif Face2 > 4 %Face = 1 or 2; V & W on face
        inc21 = (minc + 1);
        inc22 = 1;
        inc23 = ncel2;
    else %Face = 1 or 2; V & W on face
        inc21 = ncel2;
        inc22 = 1;
        inc23 = (minc + 1);
    end
    
    %Sub Model - Nodes
%     ind = ecn2(Patch2,1);
    ind = 1;
    NodeList = zeros(2,ncel2); %List of nodes that will be duplicates
    i = 1;
    for k = 1:(minc + 1)
        ind1 = (k-1)*inc21 + ind + (ind23 - 1)*inc23;
        for j = 1:(minc + 1)
            node = ind1;
            NodeList(1,i) = node; %store nodes from Face2 of Patch2 in Model2
            ind1 = ind1 + inc22;
            i = i + 1;
        end
    end
    ind1 = ind + (ind23 - 1)*inc23;
    node = ind1;
    for l = 1:3
        Corners2(1,l) = xc(l,node);
    end
    ind1 = ind1 + inc22*minc;
    node = ind1;
    for l = 1:3
        Corners2(2,l) = xc(l,node);
    end
    ind1 = minc*inc21 + ind + (ind13 - 1)*inc23;
    node = ind1;
    for l = 1:3
        Corners2(4,l) = xc(l,node);
    end
    ind1 = ind1 + inc22*minc;
    node = ind1;
    for l = 1:3
        Corners2(3,l) = xc(l,node);
    end
    
    %Check Compatibility of faces; return connectivity information and loop
    %counter increments
    [flag,inc21,inc22,start1,start2] = CheckConnectivity(Corners1,Corners2,minc,inc21,inc22);

%     if flag >= 1 %Patches are compatible; perform attachement
        %Identify Nodes in patch from Global Model that match the nodes on the
        %face of the patch in the Sub Model; these IDs must be replaced in the
        %output PatchTable
        ind = cis(1,Patch1);
        for k = 1:(minc + 1)
            ind1 = (k-1)*inc11 + ind + (ind13 - 1)*inc13;
            ind2 = (k-start1-1)*inc21 + start2*-inc22 + 1;
            for j = 1:(minc + 1)
                node = ecn(ind1);
                NodeList(2,ind2) = node;
                ind1 = ind1 + inc12;
                ind2 = ind2 + inc22;
                i = i + 1;
            end
        end
        %Add nodes from attachment face to aggregate duplicated node list
%         temp = TotalNodeList;
%         TotalNodeList = zeros(lenNL+(minc + 1)*(minc + 1),2);
%         for i = 1:lenNL
%             for j = 1:2
%                 TotalNodeList(i,j) = temp(i,j);
%             end
%         end
        for i = 1:ncel2
            for j = 1:2
                TotalNodeList(j,i+lenNL) = NodeList(j,i);
            end
        end
        %Increment number of duplicated nodes
        lenNL = lenNL + ncel2;
%     else
%         break
%     end
    
end

%Transfer data into cisout
%     cisout = zeros(numP1+numP2+1,2);
%     for i = 1:numP1
%         for j = 1:2
%             cisout(i,j) = cis1(i,j);
%         end %j
%     end %i
%Offset first indices for Sub Model by the last indices in the Global
%Model
% offset = cis(numP1+1,1) - 1;
% cis(1+numP1,1) = cis(numP1+1,1);
cis(2,1+numP1) = nel;
cis(1,2+numP1) = ncel3 + cis(1,numP1+1);
cis(2,2+numP1) = nel;

%Sort merged nodes in ascending order, remove duplicated merged nodes,
%set variables for defining new node IDs
TotalNodeList = TotalNodeList(:,1:lenNL)';
[TotalNodeList(:,1),len] = sort(TotalNodeList(:,1));
TotalNodeList(:,2) = TotalNodeList(len,2);
[TotalNodeList,lenNL] = PurgeList(TotalNodeList, lenNL, 2);
TotalNodeList = TotalNodeList';
nodestart = numnps + 1;
maxnode = TotalNodeList(1,lenNL);

%Transfer nodes from Global Model unaltered to output model
%     xs = zeros(length(xs1)+length(xc)-lenNL,4);
numnps = numnps + ncel3 - lenNL;
%     for i = 1:length(xs1)
%         for l = 1:4
%             xs(i,l) = xs1(i,l);
%         end
%     end
%Transfer from Sub Model to output model any nodes which have old IDs
%less than the least ID in NodeList
for i = 1:TotalNodeList(1,1)-1
    for l = 1:3
        xs(l,i+nodestart-1) = xc(l,i);
    end
end
%Transfer from Sub Model to output model any nodes with IDs not equal
%to any of the IDs in NodeList; performed by progressing through the
%(sorted) NodeList table and transfering nodes that fall between the
%IDs in the table
node = nodestart + TotalNodeList(1,1)-1;
j = 1;
k = TotalNodeList(1,1);
while k < maxnode
    k = k + 1;
    for i = TotalNodeList(1,j)+1:TotalNodeList(1,j+1)-1
        for l = 1:3
            xs(l,node) = xc(l,k);
        end
        node = node + 1;
        k = k + 1;
    end
    j = j + 1;
end
%Transfer from Sub Model to output model any nodes which have old IDs
%greater than the largest ID in NodeList
for i = maxnode+1:ncel3
    for l = 1:3
        xs(l,node) = xc(l,i);
    end
    node = node + 1;
end

%Adjust node references in Sub Model Patch Table

offset = nodestart - (lenNL + 1);
%The transfer of nodal coordinates was established in the above algorithm as:
%NewID = OldID + (number of nodes in Model1)
%        - (number of nodes in NodeList with IDs less than OldID)
%offset = (number of nodes in Model1 + 1) - (number of nodes in NodeList + 1)
%thus NewID = OldID + offset + (number of nodes in NodeList with IDs
%greater than OldID)

for i = 1:ncel3
    node = i;
    newnode = node + offset; %set initial new node ID
    j = lenNL; %index of highest node in NodeList
    %count number of IDs in NodeList which are greater than current node
    while j > 0 && node < TotalNodeList(1,j) 
        j = j - 1;
        newnode = newnode + 1;
    end
    if j > 0 && node == TotalNodeList(1,j)
        ecn2(i) = TotalNodeList(2,j); %assign replacement node ID
    else
        ecn2(i) = newnode; %assign new node ID
    end
end

%Transfer Patch Table data
%     ecnout = zeros(cis1(numP1+1,1)+cis2(numP2+1,1)-2,1);
node = cis(1,numP1+1)-1;
%     for i = 1:node
%         ecnout(i) = ecn1(i);
%     end
for i = 1:ncel3
    node = node + 1;
    ecn(node) = ecn2(i);
end


% Transfer data from element array to submesh array
for k = 1:cel
 cell = cell + 1;
 for j = 1:nel
    node = ixc(j,k);
    node = ecn2(node);
    ixs(j,cell) = node; %#ok<AGROW>
 end
end