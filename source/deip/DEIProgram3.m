% Program: DEIProgram3
% Tim Truster
% 08/15/2014
%
% Script to convert CG mesh (e.g. from block program output) into DG mesh
% Addition of modules to insert interfaces selectively between materials
% 
% Input: matrNodesOnElement InterTypes(nummat,nummat) for listing the interface
% material types to assign; only lower triangle is accessed; value 0 means
% interface is not inserted between those materials (including diagonal).
%
% Numbering of interface types (matI) is like:
% 1
% 2 3
% 4 5 6 ...
%
% Output:
% Actual list of used interfaces is SurfacesI and numCL, new BC arrays and
% load arrays
% Total list of interfaces is in numIfac and Interfac, although there is
% no designation on whether the interfaces are active or not. These are
% useful for assigning different DG element types to different interfaces.

% Revisions:
% 02/10/2017: adding periodic boundary conditions (PBC); uses both zipped mesh
% (collapsing linked PBC nodes into a single node and updating
% connectivity), and ghost elements (leaving original mesh, added elements
% along PBC edges that connect across the domain); both result in a 'torus'
% mesh.

% Revisions:
% 04/04/2017: adding periodic BC in the form of multi-point constraints,
% which do not require the domain to be a cube. Still uses the zipped mesh,
% but does not add ghost elements. Instead, node duplication is performed
% directly on the zipped mesh, and MPC elements are applied outside of this
% subroutine.


% flags for clearing out intermediate variables in this script; turned on
% by default
currvariables = who();
clearinter = 1;

% Test for Octave vs Matlab compatibility
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

maxel = 300;
ndm = 3;
ElementsOnNode = zeros(maxel,numnp); % Array of elements connected to nodes
ElementsOnNodeDup = zeros(maxel,numnp); % actual node ID given to element for an original node; general counterpart of NodeMat
ElementsOnNodeNum = zeros(numnp,1);
if ~exist('usePBC','var')
    usePBC = 0;
end
if usePBC % create extra arrays for zipped mesh
ElementsOnNodePBC = zeros(maxel,numnp); % Array of elements connected to nodes
ElementsOnNodePBCDup = zeros(maxel,numnp); % actual node ID given to element for an original node; general counterpart of NodeMat
ElementsOnNodePBCNum = zeros(numnp,1);
end
NodeReg = zeros(numnp,nummat);
Coordinates3 = zeros(numel*nen,ndm);
NodesOnElementCG = NodesOnElement;
if usePBC % Add space in connectivity for ghost elements
%     if nen == 4
%         NodesOnElementCG = [NodesOnElementCG zeros(numel,8-nen)
%                             zeros(numPBC*2,8)];
%         nen = 8;
%     elseif nen == 6
%         NodesOnElementCG = [NodesOnElementCG zeros(numel,8-nen)
%                             zeros(numPBC*2,8)];
%         nen = 8;
%     elseif nen == 10
%         NodesOnElementCG = [NodesOnElementCG zeros(numel,27-nen)
%                             zeros(numPBC*2,27)];
%         nen = 27;
%     elseif nen == 18
%         NodesOnElementCG = [NodesOnElementCG zeros(numel,27-nen)
%                             zeros(numPBC*2,27)];
%         nen = 27;
%     end
    NodesOnElementPBC = NodesOnElementCG; % Connectivity that is zipped, having PBC nodes condensed together
%     RegionOnElement = [RegionOnElement; zeros(numPBC*2,1)];
else
% RegionOnElement = RegionOnElement;
end
numz = 0;


%% Set up PBC data structures

if usePBC
% Lists of node pairs in convenient table format
NodesOnPBC = zeros(maxel,numnp); % nodes that are paired up
NodesOnPBCnum = zeros(numnp,1); % number of nodes with common pairing
NodesOnLink = zeros(maxel,numnp); % IDs of PBC referring to node
NodesOnLinknum = zeros(numnp,1); % number of PBC referring to node
% Handle the corners first
if usePBC == 2
    PBCList = MPCList;
    numPBC = size(PBCList,1);
else
    numPBC = size(PBCList,1);
    [TieNodes,IA] = unique(PBCList(:,3:4),'rows','first');
    if size(TieNodes,1) ~= ndm
        error('Wrong number of master pairs')
    else
        for iPBC = 1:ndm
            nodesAB = TieNodes(iPBC,1:2);
            numA = 1;
            NodesOnPBC(numA,nodesAB(1)) = nodesAB(2);
            NodesOnPBCnum(nodesAB(1)) = numA;
            numB = iPBC;
            NodesOnPBC(numB,nodesAB(2)) = nodesAB(1);
            NodesOnPBCnum(nodesAB(2)) = numB;
            numA = 1;
            NodesOnLink(numA,nodesAB(1)) = -IA(iPBC);
            NodesOnLinknum(nodesAB(1)) = numA;
            numB = iPBC;
            NodesOnLink(numB,nodesAB(2)) = -IA(iPBC);
            NodesOnLinknum(nodesAB(2)) = numB;
        end    
    end
end
% Now do the rest of the nodes
for iPBC = 1:numPBC
    nodesAB = PBCList(iPBC,1:2);
    numA = NodesOnPBCnum(nodesAB(1));
    currPBC = NodesOnPBC(1:numA,nodesAB(1));
    % make sure pair isn't already in the list for nodeA
    if isOctave
    old = ismember(nodesAB(2),currPBC);
    else
    old = ismembc(nodesAB(2),currPBC);
    end
    if ~old
        % add the pair to nodeA, expand by one
        numA = numA + 1;
        NodesOnPBC(numA,nodesAB(1)) = nodesAB(2);
        NodesOnPBCnum(nodesAB(1)) = numA;
        % Record which PBC ID refers to these nodes
        numA = NodesOnLinknum(nodesAB(1)) + 1;
        NodesOnLink(numA,nodesAB(1)) = iPBC;
        NodesOnLinknum(nodesAB(1)) = numA;
    end
    numB = NodesOnPBCnum(nodesAB(2));
    currPBC = NodesOnPBC(1:numB,nodesAB(2));
    % make sure pair isn't already in the list for nodeB
    if isOctave
    old = ismember(nodesAB(1),currPBC);
    else
    old = ismembc(nodesAB(1),currPBC);
    end
    if ~old
        % add the pair to nodeB, expand by one
        numB = numB + 1;
        NodesOnPBC(numB,nodesAB(2)) = nodesAB(1);
        NodesOnPBCnum(nodesAB(2)) = numB;
        % Record which PBC ID refers to these nodes
        numB = NodesOnLinknum(nodesAB(2)) + 1;
        NodesOnLink(numA,nodesAB(2)) = iPBC;
        NodesOnLinknum(nodesAB(2)) = numB;
    end
end
% Find additional connections for nodes with pairs of pairs; this WILL find
% the extra pairs that are usually deleted to remove linear dependence in
% the stiffness matrix, e.g. the repeated edges in 2d or 3d.
MorePBC = find(NodesOnPBCnum>1);
for i = 1:length(MorePBC)
    nodes = MorePBC(i);
    if NodesOnPBCnum(nodes) > 0 % pairs not already condensed
        notdone = 1;
        % use a tree of connectivity to get to all node pairs involving MorePBC(i)
        while notdone
            lenn = length(nodes);
            nodesS = nodes;
            for j = 1:lenn
                node = nodes(j);
                nodes = [nodes; NodesOnPBC(1:NodesOnPBCnum(node),node)];
            end
            nodesU = unique(nodes);
            if isequal(nodesS,nodesU) % starting list of node pairs matches the unique short list
                nodes = nodesU;
                notdone = 0;
            else % found some new pairs, try searching again
                nodes = nodesU;
            end
        end
        % record the maximum connected pairs into the NodesOnPBC for all
        % the connected nodes in the set; also combine together the
        % NodesOnLink at the same time
        lenn = length(nodes);
        links = [];
        for j = 1:lenn
            node = nodes(j);
            NodesOnPBCnum(node) = -(lenn - 1); % flag this node as already condensed
            NodesOnPBC(1:lenn-1,node) = setdiff(nodes,node);
            links = [links; NodesOnLink(1:NodesOnLinknum(node),node)];
        end
        links = unique(links);
        links = links(links~=0);
        lenl = length(links);
        NodesOnLink(1:lenl,nodes) = links*ones(1,lenn);
        NodesOnLinknum(nodes) = lenl;
    end
end

NodesOnPBCnum = abs(NodesOnPBCnum);

if usePBC ~= 2
    
% add in extra PBC on nodes with 3 links, so that ghost elements can always
% find their direction
MorePBC = find(NodesOnPBCnum==3);
PBCList2 = zeros(length(MorePBC)/4,4);
newPBC = numPBC;
for i = 1:length(MorePBC)
    nodes = MorePBC(i);
    if NodesOnPBCnum(nodes) > 0 % pairs not already condensed
        links = NodesOnLink(1:NodesOnLinknum(nodes),nodes);
        locPBC = PBCList(links,:);
        ind = ismember(locPBC(:,1),locPBC(:,2));
        nodesP = locPBC(ind,1:2);
        otherP = locPBC(~ind,1:2);
        if otherP(1,2) == nodesP(2)
            link = [otherP(2:-1:1,1)' locPBC(ind,3:4)];
        else
            link = [otherP(1:2,1)' locPBC(ind,3:4)];
        end
        nodes = [nodes; NodesOnPBC(1:NodesOnPBCnum(nodes),nodes)];
        newPBC = newPBC + 1;
        PBCList2(newPBC-numPBC,:) = link;
        NodesOnLink(4,nodes) = newPBC*ones(1,4);
        NodesOnPBCnum(nodes) = -4;
        NodesOnLinknum(nodes) = 4;
    end
end
NodesOnPBCnum = abs(NodesOnPBCnum);
PBCList = [PBCList; PBCList2];

else
    
% add in extra PBC on nodes with 3 links, so that ghost elements can always
% find their direction
MorePBC = find(NodesOnPBCnum>1);
PBCList2 = zeros(round(length(MorePBC)/4),2+ndm);
newPBC = numPBC;
for i = 1:length(MorePBC)
    nodes = MorePBC(i);
    numoldpairs = NodesOnPBCnum(nodes);
    if numoldpairs > 0 % pairs not already condensed
        links = NodesOnLink(1:NodesOnLinknum(nodes),nodes);
        locPBC = PBCList(links,:);
        allnodes = sort([nodes NodesOnPBC(1:numoldpairs,nodes)']);
        pairs = nchoosek(allnodes,2);
        pairs = pairs(:,2:-1:1);
        newpairs = setdiff(pairs,locPBC(:,1:2),'rows');
        numpairs = size(newpairs,1);
        for ipair = 1:numpairs
            pair = newpairs(ipair,:);
            dir_vec = Coordinates(pair(1),:)-Coordinates(pair(2),:);
            newPBC = newPBC + 1;
            PBCList2(newPBC-numPBC,:) = [pair dir_vec];
            NodesOnLink(numoldpairs+ipair,allnodes) = newPBC*ones(1,length(allnodes));
        end
        NodesOnPBCnum(allnodes) = -length(allnodes);
        NodesOnLinknum(allnodes) = (numpairs+numoldpairs);
    end
end
NodesOnPBCnum = abs(NodesOnPBCnum);
PBCList = [PBCList; PBCList2];
end

end


%% Mesh error checks

% Error checks for mesh arrays

% Nodes
if numnp < size(Coordinates,1)
    disp('Warning, number of nodes: numnp < size(Coordinates,1)')
    disp('Paused: press any key to continue')
elseif numnp > size(Coordinates,1)
    error('Error, number of nodes: numnp > size(Coordinates,1)')
end

% Elements
if numel < size(NodesOnElement,1)
    disp('Warning, number of elements: numel < size(NodesOnElement,1)')
    disp('Paused: press any key to continue')
elseif numel > size(NodesOnElement,1)
    error('Error, number of nodes: numel > size(NodesOnElement,1)')
end

if size(NodesOnElementCG,2) ~= nen
    disp('Error, connectivity: size of NodesOnElement does not match parameter nen: size(NodesOnElement,2) ~= nen')
    error('Note that FEA_Program does not use the element ID as column 1 of the array')
end

if size(Coordinates,2) ~= ndm
    disp('Warning, spatial dimensions: size of Coordinates does not match parameter ndm: size(Coordinates,2) ~= ndm')
    disp('Note that FEA_Program does not use the node ID as column 1 of the array')
    disp('Paused: press any key to continue')
end

if max(NodesOnElementCG(1:numel,1:nen)) > numnp
    error('Error, connectivity: a node in NodesOnElement exceeds the number of nodes: max(NodesOnElement(1:numel,1:nen)) > numnp')
end


%% Check that interfaces are assigned also when intraface couplers are requested
inter = find(InterTypes); % reset to logical
InterTypes(inter) = 1;
for mat = 1:nummat
    intra = InterTypes(mat,mat);
    if intra
        inter = find(InterTypes(mat+1:nummat,mat)-1);
        if ~isempty(inter)
            disp('Warning: some interfaces not requested in Intertypes and have been added')
            fprintf('intraface: %i; additional interface flags have been added\n',mat)
            InterTypes(inter+mat,mat) = 1;
        end
        inter = find(InterTypes(mat,1:mat-1)-1);
        if ~isempty(inter)
            disp('Warning: some interfaces not requested in Intertypes and have been added')
            fprintf('intraface: %i; additional interface flags have been added\n',mat)
            InterTypes(mat,inter) = 1;
        end
    end
    inter = find(InterTypes(1:mat-1,mat));
    if ~isempty(inter)
        fprintf('Warning: some entries found in the upper triangle of Intertypes(:,%i)\n',mat)
        disp('These entries are ignored for coupler insertion')
    end
end


%% Form ElementsOnNode, ElementsOnNodeNum, DG nodes and connectivity

if usePBC % also form the zipped mesh at the same time
    for elem = 1:numel

        nel = nnz(NodesOnElementPBC(elem,1:nen));
        reg = RegionOnElement(elem);

        for locN = 1:nel % Loop over local Nodes

            node = NodesOnElementPBC(elem,locN);

    %   Add element to star list, increment number of elem in star

            locE = ElementsOnNodePBCNum(node) + 1;
            ElementsOnNodePBC(locE,node) = elem;
            ElementsOnNodePBCNum(node) = locE;
            ElementsOnNodePBCDup(locE,node) = node;

              if NodesOnPBCnum(node) > 0
                nodePBC = min([node,NodesOnPBC(1,node)]); % use smallest node number as the zipped node
                NodeReg(nodePBC,reg) = nodePBC; % flag that node is in current material set
                NodesOnElementCG(elem,locN) = nodePBC;
    %   Add element to star list, increment number of elem in star
                locE = ElementsOnNodeNum(nodePBC) + 1;
                ElementsOnNode(locE,nodePBC) = elem;
                ElementsOnNodeNum(nodePBC) = locE;
                ElementsOnNodeDup(locE,nodePBC) = nodePBC;
              else
                nodePBC = node;
                NodeReg(node,reg) = node; % flag that node is in current material set
    %   Add element to star list, increment number of elem in star
                locE = ElementsOnNodeNum(nodePBC) + 1;
                ElementsOnNode(locE,nodePBC) = elem;
                ElementsOnNodeNum(nodePBC) = locE;
                ElementsOnNodeDup(locE,nodePBC) = nodePBC;
              end

        end

    end
else
    for elem = 1:numel

        nel = nnz(NodesOnElementCG(elem,1:nen));
        reg = RegionOnElement(elem);

        for locN = 1:nel % Loop over local Nodes

            node = NodesOnElementCG(elem,locN);
            NodeReg(node,reg) = node; % flag that node is in current material set

    %   Add element to star list, increment number of elem in star

            locE = ElementsOnNodeNum(node) + 1;
            ElementsOnNode(locE,node) = elem;
            ElementsOnNodeNum(node) = locE;
            ElementsOnNodeDup(locE,node) = node;


        end

    end
end

NodesOnElementDG = NodesOnElementCG;
Coordinates3(1:numnp,:) = Coordinates(1:numnp,:);

numinttype = nummat*(nummat+1)/2;
numEonF = zeros(numinttype,1);
numEonB = zeros(nummat,1);
ElementsOnFacet = zeros(8*numel,4);
FacetsOnInterface = zeros(8*numel,1);
ElementsOnBoundary = zeros(numel,2,nummat); % All exposed faces
if usePBC
FacetsOnPBC = zeros(8*numel,1); % separate list for PBC facets that are found; grouped by interface type
numEonPBC = zeros(numinttype,1);
numFPBC = 0; % record how many PBC facets are found
end
nloop = zeros(8,2);
FacetsOnElement = zeros(numel,8);
FacetsOnElementInt = zeros(numel,8);
FacetsOnNode = zeros(maxel,numnp); % for each node, which other node is connected to it by an edge; value is the nodal ID of the connecting node
FacetsOnNodeCut = FacetsOnNode; % flag for whether that edge is being cut between two nodes; 1 for cut, 0 for retain.
FacetsOnNodeInt = FacetsOnNode; % flag for materialI ID for that edge; -1 means domain edge.
FacetsOnNodeNum = zeros(numnp,1); % number of edges attached to a node.

numfac = 0;
% Templates for nodes and facets
numeT = 4;
numfnT = [3 3 3 3 3];
nloopT=[2 3 4 6 9 10
       1 3 4 7 8 10
       1 2 4 5 8 9
       1 2 3 5 6 7];
numeW = 5;
numfnW = [4 4 4 3 3];
nloopW=[2 3 5 6 8 11 14 15
       1 3 4 6 9 12 13 15
       1 2 4 5 7 10 13 14
       1 2 3 7 8 9 0 0
       4 5 6 10 11 12 0 0];
numeH = 6;
numfnH = [4 4 4 4 4 4];
nloopH=[1 4 5 8 12 16 17 20
       2 3 6 7 10 14 18 19
       1 2 5 6 9 13 17 18
       3 4 7 8 11 15 19 20
       1 2 3 4 9 10 11 12
       5 6 7 8 13 14 15 16];


%% Find all facets in mesh

for elem = 1:numel
    
    nel = nnz(NodesOnElementCG(elem,1:nen));
    
    if (nel==4||nel==10)
        
        nume = numeT;
        numfn = numfnT;
        nloop = nloopT;
        nel2 = 4;
        
    elseif (nel==6||nel==18)
        
        nume = numeW;
        numfn = numfnW;
        nloop = nloopW;
        nel2 = 6;
        
    else
        
        nume = numeH;
        numfn = numfnH;
        nloop = nloopH;
        nel2 = 8;
        
    end
    
    % Loop over edges of element
        
    for locF = 1:nume
        
        reg1 = RegionOnElement(elem); % reg1 is overwritten below, so it must be re-evaluated
        numfacnod = numfn(locF); % number of nodes on face

        nodeABCD = NodesOnElementCG(elem,nloop(locF,1:numfacnod));
        numABCD = ElementsOnNodeNum(nodeABCD);      
        ElementsOnNodeABCD = ElementsOnNode(1:max(numABCD),nodeABCD);
        
        % Determine if edge is on domain interior or boundary
        if all(numABCD>1)
            % Clean and fast way to intersect the 3 sets of elements, using
            % built-in Matlab functions; change ismembc to ismember if the
            % function is not in the standard package
            if isOctave
                twoelem = ElementsOnNodeABCD(1:numABCD(1),1);
                for j = 2:numfacnod
                    twoelem = twoelem(ismember(twoelem,ElementsOnNodeABCD(1:numABCD(j),j)));
                end
            else
                twoelem = ElementsOnNodeABCD(1:numABCD(1),1);
                for j = 2:numfacnod
                    twoelem = twoelem(ismembc(twoelem,ElementsOnNodeABCD(1:numABCD(j),j)));
                end
            end
            if length(twoelem)==2 % element interface
                facIL = 1;
                if twoelem(1)==elem
                elemA = twoelem(2);
                else
                elemA = twoelem(1);
                end
%                 if usePBC % Don't reconsider new ghost elements
%                     if elemA > numel
%                         facIL = 0;
%                     end
%                 end
            else % domain boundary
                facIL = 0;
            end
            
        else % domain boundary
            facIL = 0;
        end
        
        if facIL == 1 %element interface, add to SurfacesI
            
            % Find which slots each node occupies on elemA
            nodeAABBCCDD = zeros(4,1);
            % Other element may be a different type
            nelA = nnz(NodesOnElementCG(elemA,1:nen));
            if (nelA==4||nelA==10)
                nloopA = nloopT;
            elseif (nelA==6||nelA==18)
                nloopA = nloopW;
            else
                nloopA = nloopH;
            end
            % But the number of nodes on face won't change
            for j = 1:numfacnod
                i = 1;
                nodeAA = NodesOnElementCG(elemA,i);
                while i<nel2+1&&nodeAA~=nodeABCD(j)
                    i = i + 1;
                    nodeAA = NodesOnElementCG(elemA,i);
                end
                nodeAABBCCDD(j) = i;
            end
            
            sABCD=sort(nodeAABBCCDD(1:numfacnod))';
            for ii = 1:size(nloopA,1)
              if numfacnod == numfn(ii)
                if norm(sABCD-nloopA(ii,1:numfacnod))==0
                   locFA = ii;
                   break
                end  
              end
            end
            
            locFB = FacetsOnElement(elemA,locFA);
            
            if locFB == 0 % New edge, add to list
                
                reg2 = RegionOnElement(elemA);
                elemB = elem;
                locFB = locF;
                if reg2 > reg1 %swap the order of L and R so that L is always larger material ID
                    regA = reg1;
                    regB = reg2;
                    temp = elemB;
                    elemB = elemA;
                    elemA = temp;
                    temp = locFB;
                    locFB = locFA;
                    locFA = temp;
                else
                    regA = reg2;
                    regB = reg1;
                end
                reg1 = regA;
                reg2 = regB;
                regI = reg2*(reg2-1)/2 + reg1; % ID for material pair (row=mat2, col=mat1)
                numSI = numEonF(regI);
                numSI = numSI + 1;
                numfac = numfac + 1;
                numEonF(regI) = numSI;
                FacetsOnElement(elemB,locFB) = numfac;
                FacetsOnElement(elemA,locFA) = numfac;
                FacetsOnElementInt(elemB,locFB) = regI;
                FacetsOnElementInt(elemA,locFA) = regI;
                ElementsOnFacet(numfac,1) = elemB; %elemL
                ElementsOnFacet(numfac,2) = elemA;  %elemR
                ElementsOnFacet(numfac,3) = locFB; %facL
                ElementsOnFacet(numfac,4) = locFA;  %facR
                FacetsOnInterface(numfac) = regI;
                if usePBC
                    % Other element may be a different type
                    nelA = nnz(NodesOnElementCG(elemA,1:nen));
                    if (nelA==4||nelA==10)
                        numfnA = numfnT;
                        nloopA = nloopT;
                    elseif (nelA==6||nelA==18)
                        numfnA = numfnW;
                        nloopA = nloopW;
                    else
                        numfnA = numfnH;
                        nloopA = nloopH;
                    end
                    nodeAABBCCDD = NodesOnElementCG(elemA,nloopA(locFA,1:numfnA));
                    nodeAABBCCDD2 = NodesOnElementPBC(elemA,nloopA(locFA,1:numfnA));
                    nelB = nnz(NodesOnElementCG(elemB,1:nen));
                    if (nelB==4||nelB==10)
                        nloopB = nloopT;
                    elseif (nelB==6||nelB==18)
                        nloopB = nloopW;
                    else
                        nloopB = nloopH;
                    end
                    nodeEEFFGGHH = NodesOnElementCG(elemB,nloopB(locFB,1:numfnA));
                    nodeEEFFGGHH2 = NodesOnElementPBC(elemB,nloopB(locFB,1:numfnA));
                    
                    if any(nodeAABBCCDD~=nodeAABBCCDD2) || any(nodeEEFFGGHH~=nodeEEFFGGHH2)
                        
                      % Make sure facet is NOT an interface in the
                      % unzipped mesh, because for triangles both nodes
                      % might be in PBCList but they are on different
                      % boundary surfaces
                      numABCD = ElementsOnNodePBCNum(nodeAABBCCDD2);      
                      ElementsOnNodeABCD = ElementsOnNodePBC(1:max(numABCD),nodeAABBCCDD2);
        
                      % Determine if edge is on domain interior or boundary
                      if all(numABCD>1)
                        % Clean and fast way to intersect the 3 sets of elements, using
                        % built-in Matlab functions; change ismembc to ismember if the
                        % function is not in the standard package
                        if isOctave
                            twoelem = ElementsOnNodeABCD(1:numABCD(1),1);
                            for j = 2:numfacnod
                                twoelem = twoelem(ismember(twoelem,ElementsOnNodeABCD(1:numABCD(j),j)));
                            end
                        else
                            twoelem = ElementsOnNodeABCD(1:numABCD(1),1);
                            for j = 2:numfacnod
                                twoelem = twoelem(ismembc(twoelem,ElementsOnNodeABCD(1:numABCD(j),j)));
                            end
                        end
                        if length(twoelem)==2 % element interface
                            facIL = 1;
                            if twoelem(1)==elem
                            elemA = twoelem(2);
                            else
                            elemA = twoelem(1);
                            end
            %                 if usePBC % Don't reconsider new ghost elements
            %                     if elemA > numel
            %                         facIL = 0;
            %                     end
            %                 end
                        else % domain boundary
                            facIL = 0;
                        end

                      else % domain boundary
                        facIL = 0;
                      end
                      
                      if facIL == 0
                        cutPBC = 1;
                        numSI = numEonPBC(regI);
                        numSI = numSI + 1;
                        numEonPBC(regI) = numSI;
                        numFPBC = numFPBC + 1;
                        FacetsOnPBC(numfac) = regI;
                      else
                        cutPBC = 0;
                      end
                    else
                        cutPBC = 0;
                    end
                else
                    cutPBC = 0;
                end
                % Assign nodal edge pairs
                nodecut = InterTypes(reg2,reg1) > 0 || cutPBC;
                for j = 1:numfacnod
                    facnumABCD = FacetsOnNodeNum(nodeABCD(j)) + 1;
                    FacetsOnNodeNum(nodeABCD(j)) = facnumABCD;
                    FacetsOnNode(facnumABCD,nodeABCD(j)) = numfac;
                    FacetsOnNodeInt(facnumABCD,nodeABCD(j)) = regI;
                    if nodecut
                        FacetsOnNodeCut(facnumABCD,nodeABCD(j)) = 1;
                    end
                end
                if nel > nel2
                    % Add in the nodes along the edge too, but not the face
                    nodeABCD = NodesOnElementCG(elem,nloop(locF,numfacnod+1:2*numfacnod));
                    nodecut = InterTypes(reg2,reg1) > 0;
                    for j = 1:numfacnod
                        facnumABCD = FacetsOnNodeNum(nodeABCD(j)) + 1;
                        FacetsOnNodeNum(nodeABCD(j)) = facnumABCD;
                        FacetsOnNode(facnumABCD,nodeABCD(j)) = numfac;
                        FacetsOnNodeInt(facnumABCD,nodeABCD(j)) = regI;
                        if nodecut
                            FacetsOnNodeCut(facnumABCD,nodeABCD(j)) = 1;
                        end
                    end
                end
                
            end
            
        else %domain boundary, add to SurfaceF
            
            if usePBC
                error('domain boundaries not allowed in periodic mesh')
            end
            
            numSI = numEonB(reg1);
            numSI = numSI + 1;
            numEonB(reg1) = numSI;
            FacetsOnElement(elem,locF) = -numEonB(reg1);
            ElementsOnBoundary(numSI,1,reg1) = elem;
            ElementsOnBoundary(numSI,2,reg1) = locF;
            % Assign nodal edge pairs
            for j = 1:numfacnod
                facnumABCD = FacetsOnNodeNum(nodeABCD(j)) + 1;
                FacetsOnNodeNum(nodeABCD(j)) = facnumABCD;
                FacetsOnNode(facnumABCD,nodeABCD(j)) = numSI;
                FacetsOnNodeInt(facnumABCD,nodeABCD(j)) = -1;
            end
            if nel > nel2
                % Add in the nodes along the edge too, but not the face
                nodeABCD = NodesOnElementCG(elem,nloop(locF,numfacnod+1:2*numfacnod));
                for j = 1:numfacnod
                    facnumABCD = FacetsOnNodeNum(nodeABCD(j)) + 1;
                    FacetsOnNodeNum(nodeABCD(j)) = facnumABCD;
                    FacetsOnNode(facnumABCD,nodeABCD(j)) = numSI;
                    FacetsOnNodeInt(facnumABCD,nodeABCD(j)) = -1;
                end
            end
            
        end %facIL
        
    end %locF
    
end %elem

% Group facets according to interface type
[~,FacetsOnInterface] = sort(FacetsOnInterface(1:numfac));
FacetsOnInterfaceNum = [1; numEonF];
for int = 1:numinttype
    FacetsOnInterfaceNum(int+1) = FacetsOnInterfaceNum(int) + FacetsOnInterfaceNum(int+1);
end
numEonFPBC = numEonF;
% % The actual facet identifiers for interface type regI are then:
% locF = FacetsOnInterfaceNum(regI):(FacetsOnInterfaceNum(regI+1)-1);
% facs = FacetsOnInterface(locF);
% true = all(ElementsOnFacet(facs,1:4) == ElementsOnFacetInt(1:numEonF(regI),1:4,regI));
if usePBC % sort to find the list of PBC facets grouped by interface type
    
    [~,FacetsOnPBC] = sort(FacetsOnPBC(1:numfac));
    FacetsOnPBC = FacetsOnPBC(numfac-numFPBC+1:numfac); % delete the zeros
    FacetsOnPBCNum = [1; numEonPBC];
    for int = 1:numinttype
        FacetsOnPBCNum(int+1) = FacetsOnPBCNum(int) + FacetsOnPBCNum(int+1);
    end
    FacetsOnIntMinusPBC = zeros(numfac-numFPBC,1);
    FacetsOnIntMinusPBCNum = [1; numEonF];
    for int = 1:numinttype
        setAll = FacetsOnInterfaceNum(int):(FacetsOnInterfaceNum(int+1)-1);
        setAll = FacetsOnInterface(setAll);
        setPBC = FacetsOnPBCNum(int):(FacetsOnPBCNum(int+1)-1);
        setPBC = FacetsOnPBC(setPBC);
        setnPBC = setdiff(setAll,setPBC);
        if size(setnPBC,2) == 0
            numSI = 0;
        else
            numSI = size(setnPBC,1);
        end
        FacetsOnIntMinusPBCNum(int+1) = FacetsOnIntMinusPBCNum(int) + numSI;
        if numSI > 0
            FacetsOnIntMinusPBC(FacetsOnIntMinusPBCNum(int):(FacetsOnIntMinusPBCNum(int+1)-1)) = ...
                setnPBC;
        end
    end
    
end


%% Find material interfaces and duplicate nodes
NodesOnInterface = [];
for reg2 = 2:nummat
    for reg1 = 1:reg2-1
        
        regI = reg2*(reg2-1)/2 + reg1; % ID for material pair (row=mat2, col=mat1)
        if InterTypes(reg2,reg1) > 0
            
            % Mark nodes on the interfaces
            internodes = NodeReg(:,reg1)>0 & NodeReg(:,reg2)>0 & NodeReg(:,reg2)==NodeReg(:,reg1);
            internodes = find(internodes);
            NodesOnInterface = [NodesOnInterface; internodes];
            
        end
        
    end
end
if usePBC % Add PBC nodes into duplicating list
    NodeRegSum = sum(NodeReg,2); % from zipped nodes, ONLY the one with regions attached is in the zipped connectivity
    internodes = NodesOnPBCnum>0 & NodeRegSum>0;
    internodes = find(internodes); % add these zipped PBC nodes into the list for duplicates too
    NodesOnInterface = [NodesOnInterface; internodes];
end
NodesOnInterface = unique(NodesOnInterface);
NodesOnInterfaceNum = length(NodesOnInterface);

% Criteria: only nodes for which ALL inter-material edges involving a given
% material are being cut, then they are duplicated.
intermat2 = zeros(nummat,nummat);
for reg2 = 1:nummat
    for reg1 = 1:reg2
        intermat2(reg2,reg1) = reg2*(reg2-1)/2 + reg1;
        intermat2(reg1,reg2) = reg2*(reg2-1)/2 + reg1;
    end
end
numnp_new = numnp;
for locN = 1:NodesOnInterfaceNum
    node = NodesOnInterface(locN);
    numduplic = 0;
%     matnode = sum(NodeMat(node,:)>0);
    facnum = FacetsOnNodeNum(node);
    if facnum == 0 
        % Determine if midedge or midface node
        if nen == 10 || nen == 20 % tetrahedra, no mid-face elements
            midfaceedge = 0;
        elseif nen == 27
            elem = ElementsOnNode(1,node);
            inode = find(NodesOnElement(elem,:)==node);
            if inode > 26
                error('somehow found the center node')
            elseif inode > 20
                midfaceedge = 1;
            else
                midfaceedge = 0;
            end
        end
        if midfaceedge % midface; midedge will be handled below
          if usePBC && NodesOnPBCnum(node)>0 % handle copies of PBC nodes specially

            % Just get the old node numbers back and copy them into place
            numelems = ElementsOnNodeNum(node);
            elems = ElementsOnNode(:,node);

            % Loop over all elements attached to node
            for locE = 1:numelems

                elem = elems(locE);
                % Reset nodal ID for element in higher material ID
                nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                node_dup = NodesOnElementPBC(elem,nodeloc);
                NodesOnElementDG(elem,nodeloc) = node_dup;

            end
            
          else % regular midside node
              
            for reg2 = 2:nummat
                for reg1 = 1:reg2-1

                    regI = reg2*(reg2-1)/2 + reg1; % ID for material pair (row=mat2, col=mat1)
                    if InterTypes(reg2,reg1) > 0

                        % Duplicate nodes along material interface
                        internodes = NodeReg(node,reg1)>0 & NodeReg(node,reg2)>0 & NodeReg(node,reg2)==NodeReg(node,reg1);
                        % Namely, only find nodes that are part of materials mat1 and
                        % mat2, but ONLY if those nodal IDs have not been reset before,
                        % in which case the ID for mat1 will equal the ID for mat2.
                        % In this way, each new nodal ID is assigned to a single
                        % material.
                        if internodes

                            % Duplicate the nodal coordinates
                            numnp_new = numnp_new + 1;
                            NodeReg(node,reg2) = numnp_new;

                            % Loop over all nodes on material interface that need
                            % duplicated

                            numelems = ElementsOnNodeNum(node);
                            elems = ElementsOnNode(:,node);

                            % Loop over all elements attached to node
                            for j = 1:numelems

                                elem = elems(j);
                                reg = RegionOnElement(elem);
                                if reg == reg2
                                    % Reset nodal ID for element in higher material ID
                                    nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                                    NodesOnElementDG(elem,nodeloc) = numnp_new;
                                    Coordinates3(numnp_new,:) = Coordinates(node,:);
                                end

                            end

                        end

                    end

                end
            end
            
          end
          
        end % midface/midedge
    else % All corner nodes
        % form secs; a sec is a contiguous region of elements that is not
        % cut apart by any CZM facs
        
        % start by assuming all elements are separated
        activeSecs = ElementsOnNodeNum(node);
        Sectors = zeros(activeSecs);
        Sectors(1,1:activeSecs) = ElementsOnNode(1:activeSecs,node);
        numSecs = ones(activeSecs,1);
        
        for locF = 1:facnum
            regI = FacetsOnNodeInt(locF,node);
            if regI>0 % exclude internal facs
                
%                 intramattrue = ~isempty(find(matI==diag(intermat2),1)); I
%                 found out that the code is only putting cut in
%                 nodeedgecut for intermaterials, not the diagonal of
%                 intermat2, so no if-test is needed
                cuttrue = FacetsOnNodeCut(locF,node);
                
%                 if ~cuttrue || intramattrue % two elements should be joined into one sec, along with their neighbors currently in the sec
                if ~cuttrue % two elements should be joined into one sec, along with their neighbors currently in the sec
                    
                    facID = FacetsOnNode(locF,node);
                    elem1 = ElementsOnFacet(facID,1);
                    elem2 = ElementsOnFacet(facID,2);
                    % find the secs for each element
                    for iSec1 = 1:activeSecs
                        if numSecs(iSec1) > 0
                            elemSec = Sectors(1:numSecs(iSec1),iSec1);
                            sec1 = find(elemSec==elem1); % find is MUCH faster than ismember or ismembc
                        if sec1>0
                            break
                        end
                        end
                    end
                    for iSec2 = 1:activeSecs
                        if numSecs(iSec2) > 0
                            elemSec = Sectors(1:numSecs(iSec2),iSec2);
                            sec2 = find(elem2==elemSec);
                        if sec2>0
                            break
                        end
                        end
                    end
                    
                    % merge secs
                    if iSec2 ~= iSec1
                    Sectors(1:numSecs(iSec2)+numSecs(iSec1),iSec2) = ...
                        sort([Sectors(1:numSecs(iSec2),iSec2); Sectors(1:numSecs(iSec1),iSec1)]);
                    Sectors(:,iSec1) = 0;
                    numSecs(iSec2) = numSecs(iSec2)+numSecs(iSec1);
                    numSecs(iSec1) = 0;
                    end
                    
                end
                
            end % if external edge
            
        end
        
        % assign node IDs to each sec
        if usePBC && NodesOnPBCnum(node)>0 % handle copies of PBC nodes specially
            numcopies = NodesOnPBCnum(node)+1;
            copyonce = zeros(numcopies,1);
            for iSec = 1:activeSecs
                numelems = numSecs(iSec);
                if numelems > 0
                    
                    elems = Sectors(1:numelems,iSec);
                    elem = elems(1);
                    nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                    nodeP = NodesOnElementPBC(elem,nodeloc);
                    if nodeP > node
                        nodeloc = find(NodesOnPBC(1:numcopies-1,node)==nodeP)+1;
                    else
                        nodeloc = 1;
                    end
                    
                    if copyonce(nodeloc)
                        elems = Sectors(1:numelems,iSec);
                        numnp_new = numnp_new + 1;
                        node_dup = numnp_new;
                        Coordinates3(node_dup,:) = Coordinates(nodeP,:);

                        % Loop over all elements attached to node
                        for locE = 1:numelems

                            elem = elems(locE);
                            reg = RegionOnElement(elem);
                            NodeReg(nodeP,reg) = node_dup;
                            ind = find(elem==ElementsOnNode(1:ElementsOnNodeNum(node),node));
                            ElementsOnNodeDup(ind,node) = node_dup;
                            % Reset nodal ID for element in higher material ID
                            nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                            NodesOnElementDG(elem,nodeloc) = node_dup;

                        end
                    else
                        copyonce(nodeloc) = 1;
                        node_dup = nodeP;
                        for locE = 1:numelems

                            elem = elems(locE);
                            reg = RegionOnElement(elem);
                            NodeReg(nodeP,reg) = node_dup;
                            ind = find(elem==ElementsOnNode(1:ElementsOnNodeNum(node),node));
                            ElementsOnNodeDup(ind,node) = node_dup;
                            % Reset nodal ID for element in higher material ID
                            nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                            NodesOnElementDG(elem,nodeloc) = node_dup;

                        end
                    end
                end
            end
        else % regular node
            copyonce = 0;
            for iSec = 1:activeSecs
                numelems = numSecs(iSec);
                if numelems > 0
                    if copyonce
                        elems = Sectors(1:numelems,iSec);
                        numnp_new = numnp_new + 1;
                        Coordinates3(numnp_new,:) = Coordinates(node,:);

                        % Loop over all elements attached to node
                        for j = 1:numelems

                            elem = elems(j);
                            reg = RegionOnElement(elem);
                            NodeReg(node,reg) = numnp_new;
                            ind = find(elem==ElementsOnNode(1:ElementsOnNodeNum(node),node));
                            ElementsOnNodeDup(ind,node) = numnp_new;
                            % Reset nodal ID for element in higher material ID
                            nodeloc = find(NodesOnElementCG(elem,1:nen)==node);
                            NodesOnElementDG(elem,nodeloc) = numnp_new;

                        end
                    else
                        copyonce  = 1;
                    end
                end
            end
        end
        
    end %if facnum
    
end


%% Form interfaces between all edges on interiors of specified material regions
if norm(diag(InterTypes)) > 0
    
NodesOnElementCG = NodesOnElementDG; % Update connectivities to reflect presence of interfaces
numnp_I = numnp_new;
ElementsOnNodeNum2 = zeros(numnp_I,1);
NodeCGDG = zeros(maxel,numnp_I); % Mapping of DG nodal IDs to CG nodal IDs
NodeCGDG(1,1:numnp_I) = 1:numnp_I;

for reg2 = 1:nummat
    
    regI = reg2*(reg2+1)/2; % ID for material pair (row=mat2, col=mat1)
    if InterTypes(reg2,reg2) > 0
        
        % Find elements belonging to material mat2
        elems = RegionOnElement(:)==reg2;
        elems = find(elems);
        numelems = length(elems);
        
        % Loop over elements in material mat2, explode all the elements
        for i = 1:numelems

            elem = elems(i);
            nel = nnz(NodesOnElementDG(elem,1:nen));

            for locN = 1:nel % Loop over local Nodes

                node = NodesOnElementCG(elem,locN);
                if ElementsOnNodeNum2(node) == 0 % only duplicate nodes after the first time encountered
                    ElementsOnNodeNum2(node) = 1;
                    NodeCGDG(1,node) = node;    %node in DG mesh with same coordinates
                else
                    numnp_new = numnp_new + 1;
                    NodesOnElementDG(elem,locN) = numnp_new;
                    Coordinates3(numnp_new,:) = Coordinates3(node,:);
                    oldnode = NodesOnElement(elem,locN);
                    ind = find(elem==ElementsOnNode(1:ElementsOnNodeNum(oldnode),oldnode));
                    ElementsOnNodeDup(ind,oldnode) = numnp_new;

                    %   Add element to star list, increment number of elem in star

                    locE = ElementsOnNodeNum2(node) + 1;
                    ElementsOnNodeNum2(node) = locE;
                    NodeCGDG(locE,node) = numnp_new;    %node in DG mesh with same coordinates
                end

            end
            
        end
        
    end
    
end

inds = find(ElementsOnNodeNum2==0);
ElementsOnNodeNum2(inds) = 1;

else % no DG in interior of any materials
    
numnp_I = numnp_new;
ElementsOnNodeNum2 = ones(numnp_I,1); % set one copy of duplicated nodes per material ID so that BC expanding subroutine works properly.
NodeCGDG = zeros(maxel,numnp_I); % Mapping of DG nodal IDs to CG nodal IDs
NodeCGDG(1,1:numnp_I) = 1:numnp_I;

end
    

%% Form final lists of DG interfaces
numnp = numnp_new;
Coordinates = Coordinates3(1:numnp,:);
numCL = 0;
NodesOnElement = NodesOnElementDG;
if usePBC
NodesOnElementCG = NodesOnElementPBC;
end
if usePBC == 2
    MPCList = PBCList;
end


%% Clear out intermediate variables
if clearinter
KeyList = {'numEonB','numEonF','ElementsOnBoundary','numSI','ElementsOnFacet',...
'ElementsOnNode','ElementsOnNodeA','ElementsOnNodeB','ElementsOnNodeDup',...
'ElementsOnNodeNum','numfac','ElementsOnNodeNum2','numinttype','FacetsOnElement',...
'FacetsOnElementInt','FacetsOnInterface','FacetsOnInterfaceNum','FacetsOnNode',...
'FacetsOnNodeCut','FacetsOnNodeInt','FacetsOnNodeNum','NodeCGDG','NodeReg',...
'NodesOnElementCG','NodesOnElementDG','NodesOnInterface','NodesOnInterfaceNum',...
'numCL','maxel','numelPBC','RegionOnElementDG','numPBC','TieNodesNew','TieNodes',...
'NodesOnPBC','NodesOnPBCnum','NodesOnLink','NodesOnLinknum','numEonPBC',...
'FacetsOnPBC','FacetsOnPBCNum','FacetsOnIntMinusPBC','FacetsOnIntMinusPBCNum','MPCList'};
if isOctave
    wholething = [{'-x'},KeyList,currvariables'];
    clear(wholething{:})
else
    wholething = [{'-except'},KeyList,currvariables'];
    clearvars(wholething{:})
end
end
