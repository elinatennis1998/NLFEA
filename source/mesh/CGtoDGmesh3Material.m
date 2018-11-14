% Tim Truster
% 08/15/2014
%
% Script to convert CG mesh (e.g. from block program output) into DG mesh
% Addition of modules to insert interfaces selectively between materials
% 
% Input: matrix InterTypes(nummat,nummat) for listing the interface
% material types to assign; only lower triangle is accessed; value 0 means
% interface is not inserted between those materials (including diagonal).
%
% Numbering of interface types (matI) is like:
% 1
% 2 3
% 4 5 6 ...

maxel = 300;12;
ndm = 3;
epatch = zeros(maxel,numnp); % Array of elements connected to nodes
epnum = zeros(numnp,1);
NodeMat = zeros(numnp,nummat);
Coordinates3 = zeros(numel*nen,ndm);
ixCG = ix;
ixDG = ixCG;
numz = 0;


%% Form epatch, epnum, DG nodes and connectivity

for elem = 1:numel

    nel = nnz(ixCG(elem,1:nen));
    mat = ix(elem,nen1);
    ixDG(elem,nen1) = mat;
    
    for k = 1:nel % Loop over local Nodes

%         numz = numz + 1;
%         ixDG(elem,k) = numz;
        node = ixCG(elem,k);
        NodeMat(node,mat) = node; % flag that node is in current material set
%         Coordinates3(numz,:) = Coordinates(node,:);

%   Add element to star list, increment number of elem in star

        q = epnum(node) + 1;
        epatch(q,node) = elem;		%epatch(nel,numnp)
        epnum(node) = q;			%epnum(numnp)
%         NodeCGDG(q,node) = numz;    %node in DG mesh with same coordinates


    end %k
    
end %j

Coordinates3(1:numnp,:) = Coordinates(1:numnp,:);

numinttype = nummat*(nummat+1)/2;
numIEdge = zeros(numinttype,1);
numDEdge = zeros(nummat,1);
numComp = 0;
InterEdge = zeros(numel*4,8,numinttype);
DomEdge = zeros(numel,2,nummat); % All exposed faces
nloop = zeros(8,2);
elemedge = zeros(numel,8);
elemedgemat = zeros(numel,8);
nodeedge = zeros(maxel,numnp); % for each node, which other node is connected to it by an edge; value is the nodal ID of the connecting node
nodeedgecut = nodeedge; % flag for whether that edge is being cut between two nodes; 1 for cut, 0 for retain.
nodeedgemat = nodeedge; % flag for materialI ID for that edge; -1 means domain edge.
nodeedgenum = zeros(numnp,1); % number of edges attached to a node.


%% Find all edges in mesh

for elem = 1:numel
    
    nel = nnz(ixCG(elem,1:nen));
    matL = ixCG(elem,nen1);
    
    if (nel==4||nel==10)
        
        nume = 4;
        nloop=[2 3 4 6 9 10
               1 3 4 7 8 10
               1 2 4 5 8 9
               1 2 3 5 6 7];
        
    else
        nume = 6;
        nloop=[1 4 5 8 12 16 17 20
               2 3 6 7 10 14 18 19
               1 2 5 6 9 13 17 18
               3 4 7 8 11 15 19 20
               1 2 3 4 9 10 11 12
               5 6 7 8 13 14 15 16];
        
    end
    
    %     Loop over edges of element
    if (nel==4||nel==10)
        
    for edge = 1:nume

        nodeA = ixCG(elem,nloop(edge,1));
        nodeB = ixCG(elem,nloop(edge,2));
        nodeC = ixCG(elem,nloop(edge,3));
        numA = epnum(nodeA); %epnum check if the node is in the interface
        numB = epnum(nodeB);
        numC = epnum(nodeC);     
        epatchA = epatch(:,nodeA);
        epatchB = epatch(:,nodeB);
        epatchC = epatch(:,nodeC);
        
        % Determine if edge is on domain interior or boundary
        if (numA>1&&numB>1&&numC>1)

            i = 0;
            notdone = 1;
            
            % Try to find another element that has both nodeA and nodeB
            while notdone==1
                
                i = i + 1;
                elemA = epatchA(i);
                elemB = 0;
                elemC = 0;
                if elemA == elem 
                    j = numB + 1;
                    k = numC + 1;
                else
                    j = 0;
                    k = 0;
                end
                
                while j<numB&&elemB~=elemA
                    
                    j = j + 1;
                    elemB = epatchB(j);
                    
                end
                
                while k<numC&&elemC~=elemA
                    
                    k = k + 1;
                    elemC = epatchC(k);
                    
                end
                
                if elemA==elemB&&elemA==elemC||i==numA
                    notdone = 0;
                end
                
            end
            
            if elemA==elemB&&elemA==elemB&&elemA==elemC % element interface
                edgeIL = 1;
            else % domain boundary
                edgeIL = 0;
            end
            
        else % domain boundary
            edgeIL = 0;
        end
        
        if edgeIL == 1 %element interface, add to SurfacesI
            
            % Find which slots nodeA and nodeB occupy on elemA
            i = 1;
            nodeAA = ixCG(elemA,i);
            while i<10&&nodeAA~=nodeA
                i = i + 1;
                nodeAA = ixCG(elemA,i);
            end
            nodeAA = i;
            i = 1;
            nodeBB = ixCG(elemA,i);
            while i<10&&nodeBB~=nodeB
                i = i + 1;
                nodeBB = ixCG(elemA,i);
            end
            nodeBB = i;
            i = 1;
            nodeCC = ixCG(elemA,i);
            while i<10&&nodeCC~=nodeC
                i = i + 1;
                nodeCC = ixCG(elemA,i);
            end
            nodeCC = i;
            i = 1;
            sABC=sort([nodeAA nodeBB nodeCC]);
           for ii = 1:size(nloop,1)
               if norm(sABC-nloop(ii,1:3))==0
                  edgeA = ii;
               end
           end
            
            edgeB = elemedge(elemA,edgeA);
            
            if edgeB == 0 % New edge, add to list
                
                matR = ixCG(elemA,nen1);
                elemB = elem;
                edgeB = edge;
                if matR > matL %swap the order of L and R so that L is always larger material ID
                    matA = matL;
                    matB = matR;
                    mat1 = elemB;
                    elemB = elemA;
                    elemA = mat1;
                    mat1 = edgeB;
                    edgeB = edgeA;
                    edgeA = mat1;
                else
                    matA = matR;
                    matB = matL;
                end
                mat1 = matA;
                mat2 = matB;
                matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
                numSI = numIEdge(matI);
                numSI = numSI + 1;
                numIEdge(matI) = numSI;
                elemedge(elemB,edgeB) = numSI;
                elemedge(elemA,edgeA) = numSI;
                elemedgemat(elemB,edgeB) = matI;
                elemedgemat(elemA,edgeA) = matI;
                InterEdge(numSI,1,matI) = 0;
                InterEdge(numSI,2,matI) = 0;
                InterEdge(numSI,3,matI) = 0;
                InterEdge(numSI,4,matI) = 0;
                InterEdge(numSI,5,matI) = elemB; %elemL
                InterEdge(numSI,6,matI) = elemA;  %elemR
                InterEdge(numSI,7,matI) = edgeB; %edgeL
                InterEdge(numSI,8,matI) = edgeA;  %edgeR
                % Assign nodal edge pairs
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = matI;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = matI;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = matI;
                if nen > 4
                % Add in the nodes along the edge too, but not the face
                nodeA = ixCG(elem,nloop(edge,4));
                nodeB = ixCG(elem,nloop(edge,5));
                nodeC = ixCG(elem,nloop(edge,6));
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = matI;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = matI;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = matI;
                end
                
            end
            
        else %domain boundary, add to SurfaceF
            
            numSI = numDEdge(matL);
            numSI = numSI + 1;
            numDEdge(matL) = numSI;
            elemedge(elem,edge) = numDEdge(matL);
            DomEdge(numSI,1,matL) = ixCG(elem,nloop(edge,2));
            DomEdge(numSI,2,matL) = ixCG(elem,nloop(edge,1));
            DomEdge(numSI,3,matL) = elem;
            DomEdge(numSI,4,matL) = edge;
            % Assign nodal edge pairs
            edgnum = nodeedgenum(nodeA) + 1;
            nodeedgenum(nodeA) = edgnum;
%             nodeedge(edgnum,nodeA) = nodeB;
            nodeedgemat(edgnum,nodeA) = -1;
            edgnum = nodeedgenum(nodeB) + 1;
            nodeedgenum(nodeB) = edgnum;
%             nodeedge(edgnum,nodeB) = nodeA;
            nodeedgemat(edgnum,nodeB) = -1;
            edgnum = nodeedgenum(nodeC) + 1;
            nodeedgenum(nodeC) = edgnum;
%             nodeedge(edgnum,nodeB) = nodeA;
            nodeedgemat(edgnum,nodeC) = -1;
            if nen > 4
                % Add in the nodes along the edge too, but not the face
                nodeA = ixCG(elem,nloop(edge,4));
                nodeB = ixCG(elem,nloop(edge,5));
                nodeC = ixCG(elem,nloop(edge,6));
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = -1;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = -1;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = -1;
            end
            
        end
    
    end
        
    else
        
    for edge = 1:nume

        nodeA = ixCG(elem,nloop(edge,1));
        nodeB = ixCG(elem,nloop(edge,2));
        nodeC = ixCG(elem,nloop(edge,3));
        nodeD = ixCG(elem,nloop(edge,4));
        numA = epnum(nodeA); %epnum check if the node is in the interface
        numB = epnum(nodeB);
        numC = epnum(nodeC);
        numD = epnum(nodeD);        
        epatchA = epatch(:,nodeA);
        epatchB = epatch(:,nodeB);
        epatchC = epatch(:,nodeC);
        epatchD = epatch(:,nodeD);
        
        % Determine if edge is on domain interior or boundary
        if (numA>1&&numB>1&&numC>1&&numD>1)

            i = 0;
            notdone = 1;
            
            % Try to find another element that has both nodeA and nodeB
            while notdone==1
                
                i = i + 1;
                elemA = epatchA(i);
                elemB = 0;
                elemC = 0;
                elemD = 0;
                if elemA == elem 
                    j = numB + 1;
                    k = numC + 1;
                    l = numD + 1;
                else
                    j = 0;
                    k = 0;
                    l = 0;
                end
                
                while j<numB&&elemB~=elemA
                    
                    j = j + 1;
                    elemB = epatchB(j);
                    
                end
                
                while k<numC&&elemC~=elemA
                    
                    k = k + 1;
                    elemC = epatchC(k);
                    
                end                

                while l<numD&&elemD~=elemA
                    
                    l = l + 1;
                    elemD = epatchD(l);
                    
                end   
                
                if elemA==elemB&&elemA==elemC&&elemA==elemD||i==numA
                    notdone = 0;
                end
                
            end
            
            if elemA==elemB&&elemA==elemB&&elemA==elemC&&elemA==elemD % element interface
                edgeIL = 1;
            else % domain boundary
                edgeIL = 0;
            end
            
        else % domain boundary
            edgeIL = 0;
        end
        
        if edgeIL == 1 %element interface, add to SurfacesI
            
            % Find which slots nodeA and nodeB occupy on elemA
            i = 1;
            nodeAA = ixCG(elemA,i);
            while i<9&&nodeAA~=nodeA
                i = i + 1;
                nodeAA = ixCG(elemA,i);
            end
            nodeAA = i;
            i = 1;
            nodeBB = ixCG(elemA,i);
            while i<9&&nodeBB~=nodeB
                i = i + 1;
                nodeBB = ixCG(elemA,i);
            end
            nodeBB = i;
            i = 1;
            nodeCC = ixCG(elemA,i);
            while i<9&&nodeCC~=nodeC
                i = i + 1;
                nodeCC = ixCG(elemA,i);
            end
            nodeCC = i;
            i = 1;
            nodeDD = ixCG(elemA,i);
            while i<9&&nodeDD~=nodeD
                i = i + 1;
                nodeDD = ixCG(elemA,i);
            end
            nodeDD = i;
            
            sABCD=sort([nodeAA nodeBB nodeCC nodeDD]);
            for ii = 1:size(nloop,1)
                if norm(sABCD-nloop(ii,1:4))==0
                   edgeA = ii;
                end     
            end
            
            edgeB = elemedge(elemA,edgeA);
            
            if edgeB == 0 % New edge, add to list
                
                matR = ixCG(elemA,nen1);
                elemB = elem;
                edgeB = edge;
                if matR > matL %swap the order of L and R so that L is always larger material ID
                    mat1 = matL;
                    matL = matR;
                    matR = mat1;
                    mat1 = elemB;
                    elemB = elemA;
                    elemA = mat1;
                    mat1 = edgeB;
                    edgeB = edgeA;
                    edgeA = mat1;
                end
                mat1 = matR;
                mat2 = matL;
                matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
                numSI = numIEdge(matI);
                numSI = numSI + 1;
                numIEdge(matI) = numSI;
                elemedge(elemB,edgeB) = numSI;
                elemedge(elemA,edgeA) = numSI;
                elemedgemat(elemB,edgeB) = matI;
                elemedgemat(elemA,edgeA) = matI;
                InterEdge(numSI,1,matI) = 0;
                InterEdge(numSI,2,matI) = 0;
                InterEdge(numSI,3,matI) = 0;
                InterEdge(numSI,4,matI) = 0;
                InterEdge(numSI,5,matI) = elemB; %elemL
                InterEdge(numSI,6,matI) = elemA;  %elemR
                InterEdge(numSI,7,matI) = edgeB; %edgeL
                InterEdge(numSI,8,matI) = edgeA;  %edgeR
                % Assign nodal edge pairs
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = matI;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = matI;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = matI;
                edgnum = nodeedgenum(nodeD) + 1;
                nodeedgenum(nodeD) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeD) = matI;
                if nen > 8
                % Add in the nodes along the edge too, but not the face
                nodeA = ixCG(elem,nloop(edge,5));
                nodeB = ixCG(elem,nloop(edge,6));
                nodeC = ixCG(elem,nloop(edge,7));
                nodeD = ixCG(elem,nloop(edge,8));
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = matI;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = matI;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = matI;
                edgnum = nodeedgenum(nodeD) + 1;
                nodeedgenum(nodeD) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeD) = matI;
                end
                
            end
            
        else %domain boundary, add to SurfaceF
            
            numSI = numDEdge(matL);
            numSI = numSI + 1;
            numDEdge(matL) = numSI;
            elemedge(elem,edge) = numDEdge(matL);
            DomEdge(numSI,1,matL) = ixCG(elem,nloop(edge,2));
            DomEdge(numSI,2,matL) = ixCG(elem,nloop(edge,1));
            DomEdge(numSI,3,matL) = elem;
            DomEdge(numSI,4,matL) = edge;
                % Assign nodal edge pairs
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = -1;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = -1;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = -1;
                edgnum = nodeedgenum(nodeD) + 1;
                nodeedgenum(nodeD) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeD) = -1;
                if nen > 8
                % Add in the nodes along the edge too, but not the face
                nodeA = ixCG(elem,nloop(edge,5));
                nodeB = ixCG(elem,nloop(edge,6));
                nodeC = ixCG(elem,nloop(edge,7));
                nodeD = ixCG(elem,nloop(edge,8));
                edgnum = nodeedgenum(nodeA) + 1;
                nodeedgenum(nodeA) = edgnum;
%                 nodeedge(edgnum,nodeA) = nodeB;
                nodeedgemat(edgnum,nodeA) = -1;
                edgnum = nodeedgenum(nodeB) + 1;
                nodeedgenum(nodeB) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeB) = -1;
                edgnum = nodeedgenum(nodeC) + 1;
                nodeedgenum(nodeC) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeC) = -1;
                edgnum = nodeedgenum(nodeD) + 1;
                nodeedgenum(nodeD) = edgnum;
%                 nodeedge(edgnum,nodeB) = nodeA;
                nodeedgemat(edgnum,nodeD) = -1;
                end
            
        end
    
    end
    
    end
    
end


%% Find material interfaces and duplicate nodes
internodelist = [];
for mat2 = 2:nummat
    for mat1 = 1:mat2-1
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(mat2,mat1) > 0
            
            % Mark edges as being cut
            nodeedgecut(nodeedgemat==matI) = 1;
            % Mark nodes on the interfaces
            internodes = NodeMat(:,mat1)>0 & NodeMat(:,mat2)>0 & NodeMat(:,mat2)==NodeMat(:,mat1);
            internodes = find(internodes);
            internodelist = [internodelist; internodes];
            
        end
        
    end
end
internodelist = unique(internodelist);
numinternodes = length(internodelist);
% Now actually duplicate the nodes
% Criteria: only nodes for which ALL inter-material edges involving a given
% material are being cut, then they are duplicated.
intermat2 = zeros(nummat,nummat);
for mat2 = 1:nummat
    for mat1 = 1:mat2
        intermat2(mat2,mat1) = mat2*(mat2-1)/2 + mat1;
        intermat2(mat1,mat2) = mat2*(mat2-1)/2 + mat1;
    end
end
numnp_new = numnp;
for k = 1:numinternodes
    node = internodelist(k);
    numduplic = 0;
    matnode = sum(NodeMat(node,:)>0);
    edgnum = nodeedgenum(node);
    if edgnum == 0 
        % Determine if midedge or midface node
        if nen == 10 || nen == 20 % tetrahedra, no mid-face elements
            midfaceedge = 0;
        elseif nen == 27
            elem = epatch(1,node);
            inode = find(ix(elem,:)==node);
            if inode > 26
                error('somehow found the center node')
            elseif inode > 20
                midfaceedge = 1;
            else
                midfaceedge = 0;
            end
        end
        if midfaceedge % midface; midedge will be handled below
        for mat2 = 2:nummat
            for mat1 = 1:mat2-1
                
                matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
                if InterTypes(mat2,mat1) > 0
                    
                    % Duplicate nodes along material interface
                    internodes = NodeMat(node,mat1)>0 & NodeMat(node,mat2)>0 & NodeMat(node,mat2)==NodeMat(node,mat1);
                    % Namely, only find nodes that are part of materials mat1 and
                    % mat2, but ONLY if those nodal IDs have not been reset before,
                    % in which case the ID for mat1 will equal the ID for mat2.
                    % In this way, each new nodal ID is assigned to a single
                    % material.
                    if internodes
                    
                        % Duplicate the nodal coordinates
                        numnp_new = numnp_new + 1;
                        Coordinates3(numnp_new,:) = Coordinates(node,:);
                        NodeMat(node,mat2) = numnp_new;

                        % Loop over all nodes on material interface that need
                        % duplicated
                                
                        numelems = epnum(node);
                        elems = epatch(:,node);
                        
                        % Loop over all elements attached to node
                        for j = 1:numelems
                            
                            elem = elems(j);
                            mat = ixCG(elem,nen1);
                            if mat == mat2
                                % Reset nodal ID for element in higher material ID
                                nodeloc = find(ixCG(elem,1:nen)==node);
                                ixDG(elem,nodeloc) = numnp_new;
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
        end
        end % midface/midedge
    else % All corner nodes
    for ma = 1:nummat
        
        if NodeMat(node,ma) > 0
        
        % Make list of all inter-material types that involve that material
        intermat = intermat2(:,ma);
        intermat = intermat([1:ma-1 ma+1:nummat]); %remove the intramaterial one
        
        % Check if all such edges of those types are being cut
        cuttrue = 1;
        edge = 0;
        while cuttrue && edge < edgnum
            edge = edge + 1;
            matI = nodeedgemat(edge,node);
            if ismember(matI,intermat)
                cuttrue = nodeedgecut(edge,node);
            end
        end
        
        % If true, then give a new number, but only if numduplic < matnode
        % - 1
        if cuttrue && numduplic < matnode - 1
            numduplic = numduplic + 1;
            numnp_new = numnp_new + 1;
            Coordinates3(numnp_new,:) = Coordinates(node,:);
            NodeMat(node,ma) = numnp_new;
            numelems = epnum(node);
            elems = epatch(:,node);

            % Loop over all elements attached to node
            for j = 1:numelems

                elem = elems(j);
                mat = ixCG(elem,nen1);
                if mat == ma
                    % Reset nodal ID for element in higher material ID
                    nodeloc = find(ixCG(elem,1:nen)==node);
                    ixDG(elem,nodeloc) = numnp_new;
                end

            end
        end
        
        end
        
    end %for ma
    end %if edgnum
end


%% Form interfaces between all edges on interiors of specified material regions
if norm(diag(InterTypes)) > 0
    
ixCG = ixDG; % Update connectivities to reflect presence of interfaces
numnp_I = numnp_new;
epnum2 = zeros(numnp_I,1);
NodeCGDG = zeros(maxel,numnp_I); % Mapping of DG nodal IDs to CG nodal IDs
NodeCGDG(1,1:numnp_I) = 1:numnp_I;

for mat2 = 1:nummat
    
    matI = mat2*(mat2+1)/2; % ID for material pair (row=mat2, col=mat1)
    if InterTypes(mat2,mat2) > 0
        
        % Find elements belonging to material mat2
        elems = ixDG(:,nen1)==mat2;
        elems = find(elems);
        numelems = length(elems);
        
        % Loop over elements in material mat2, explode all the elements
        for i = 1:numelems

            elem = elems(i);
            nel = nnz(ixDG(elem,1:nen));

            for k = 1:nel % Loop over local Nodes

                node = ixCG(elem,k);
                if epnum2(node) == 0 % only duplicate nodes after the first time encountered
                    epnum2(node) = 1;
                    NodeCGDG(1,node) = node;    %node in DG mesh with same coordinates
                else
                    numnp_new = numnp_new + 1;
                    ixDG(elem,k) = numnp_new;
                    Coordinates3(numnp_new,:) = Coordinates3(node,:);

                    %   Add element to star list, increment number of elem in star

                    q = epnum2(node) + 1;
                    epnum2(node) = q;
                    NodeCGDG(q,node) = numnp_new;    %node in DG mesh with same coordinates
                end

            end %k
            
        end
        
    end
    
end

inds = find(epnum2==0);
epnum2(inds) = 1;

else % no DG in interior of any materials
    
numnp_I = numnp_new;
epnum2 = ones(numnp_I,1); % set one copy of duplicated nodes per material ID so that BC expanding subroutine works properly.
NodeCGDG = zeros(maxel,numnp_I); % Mapping of DG nodal IDs to CG nodal IDs
NodeCGDG(1,1:numnp_I) = 1:numnp_I;

end
    

%% Form final lists of DG interfaces
numnp = numnp_new;
Coordinates = Coordinates3(1:numnp,:);
numCL = 0;
numFL = 0;
SurfacesI = [];
for mat2 = 1:nummat
    for mat1 = 1:nummat
        
        matI = mat2*(mat2-1)/2 + mat1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(mat2,mat1) > 0
            numCL = numCL + numIEdge(matI);
            SurfacesI = [SurfacesI squeeze(InterEdge(1:numIEdge(matI),1:8,matI))']; %#ok<AGROW>
        end
    end
end
SurfacesI = SurfacesI';

clear('numD') % remove so that NL_FEA_Program doesn't find it
% SurfacesI = SurfacesI(1:numCL,:);
% SurfacesF = SurfacesF(1:numFL,:);

%% Update boundary conditions
if exist('numBC','var') && numBC > 0
    NodeBCCG = NodeBC;
    numBCCG = numBC;
    [NodeBC,numBC] = updateCGDGnodelist2(maxel,epnum2,NodeCGDG,NodeMat,NodeBCCG,numBCCG);
end
if exist('numBCnp','var') && numBCnp > 0
    NodeBCnpCG = NodeBCnp;
    numBCnpCG = numBCnp;
    [NodeBCnp,numBCnp] = updateCGDGnodelist2(maxel,epnum2,NodeCGDG,NodeMat,NodeBCnpCG,numBCnpCG);
end
% Update loads
% NODAL LOADS PROBABLY SHOULDN'T BE DUPLICATED!!!!
if exist('numNodalF','var') && numNodalF > 0
    NodeLoadCG = NodeLoad;
    numNodalFCG = numNodalF;
    [NodeLoad,numNodalF] = updateCGDGnodelist2(maxel,epnum2,NodeCGDG,NodeMat,NodeLoadCG,numNodalFCG);
end
if exist('numNodalFnp','var') && numNodalFnp > 0
    NodeLoadnpCG = NodeLoadnp;
    numNodalFnpCG = numNodalFnp;
    [NodeLoadnp,numNodalFnp] = updateCGDGnodelist2(maxel,epnum2,NodeCGDG,NodeMat,NodeLoadnpCG,numNodalFnpCG);
end

% For 3d, surface loads don't use node IDs, so we don't need to update them
if exist('numSL','var') && numSL > 0
    SurfacesLCG = SurfacesL;
    numSLCG = numSL;
    [SurfacesL,numSL] = updateCGDGelemlist(ndm,nen,ixCG,ixDG,SurfacesLCG,numSLCG);
end
if exist('numSLnp','var') && numSLnp > 0
    SurfacesLnpCG = SurfacesLnp;
    numSLCGnp = numSLnp;
    [SurfacesLnp,numSLnp] = updateCGDGelemlist(ndm,nen,ixCG,ixDG,SurfacesLnpCG,numSLCGnp);
end
if exist('numBD','var') && numBD > 0
    SurfacesDCG = SurfacesD;
    numBDCG = numBD;
    [SurfacesD,numBD] = updateCGDGelemlist(ndm,nen,ixCG,ixDG,SurfacesDCG,numBDCG);
end
if exist('numUSL','var') && numUSL > 0
    SurfacesUCG = SurfacesU;
    numUSLCG = numUSL;
    [SurfacesU,numUSL] = updateCGDGelemlist(ndm,nen,ixCG,ixDG,SurfacesUCG,numUSLCG);
end

% Final connectivity
ix = ixDG;
