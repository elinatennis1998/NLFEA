% Tim Truster
% 06/15/2011
%
% Script to convert CG mesh (e.g. from block program output) into DG mesh

maxel = 12;
ndm = 2;
epatch = zeros(maxel,numnp); % Array of elements connected to nodes
epnum = zeros(numnp,1);
NodeCGDG = epatch; % Array of CG nodes and the corresponding DG nodes
Coordinates3 = zeros(numel*nen,ndm);
NodesOnElementDG = zeros(numel,nen);
RegionOnElementDG = RegionOnElement;
NodesOnElementCG = NodesOnElement;
RegionOnElementCG = RegionOnElement;
numz = 0;

% Form epatch, epnum, DG nodes and connectivity

for elem = 1:numel

    nel = nnz(NodesOnElementCG(elem,1:nen));
    
    for k = 1:nel % Loop over local Nodes

        numz = numz + 1;
        NodesOnElementDG(elem,k) = numz;
        node = NodesOnElementCG(elem,k);
        Coordinates3(numz,:) = Coordinates(node,:);

%   Add element to star list, increment number of elem in star

        q = epnum(node) + 1;
        epatch(q,node) = elem;		%epatch(nel,numnp)
        epnum(node) = q;			%epnum(numnp)
        NodeCGDG(q,node) = numz;    %node in DG mesh with same coordinates


    end %k
    
end %j

NodesOnElement = NodesOnElementDG;
numnp = numz;
Coordinates = Coordinates3(1:numnp,:);

numCL = 0;
numFL = 0;
numComp = 0;
SurfacesI = zeros(numel*4,8);
SurfacesF = zeros(numel,7); % All exposed faces
nloop = zeros(4,2);
elemedge = zeros(numel,4);

% Form SurfaceI, SurfaceL

for elem = 1:numel
    
    nel = nnz(NodesOnElementCG(elem,1:nen));
    
    if (nel==3||nel==6)
        
        nume = 3;
        for i=1:3
            nloop(i,1) = i;
            nloop(i,2) = i+1;
        end
        nloop(3,2) = 1;
        
    else
        
        nume = 4;
        for i=1:4
            nloop(i,1) = i;
            nloop(i,2) = i+1;
        end
        nloop(4,2) = 1;
        
    end
    
    %     Loop over edges of element
    for edge = 1:nume

        nodeA = NodesOnElementCG(elem,nloop(edge,1));
        nodeB = NodesOnElementCG(elem,nloop(edge,2));
        numA = epnum(nodeA);
        numB = epnum(nodeB);
        epatchA = epatch(:,nodeA);
        epatchB = epatch(:,nodeB);
        
        % Determine if edge is on domain interior or boundary
        if (numA>1&&numB>1)
            
            i = 0;
            notdone = 1;
            
            % Try to find another element that has both nodeA and nodeB
            while notdone==1
                
                i = i + 1;
                elemA = epatchA(i);
                elemB = 0;
                if elemA == elem
                    j = numB + 1;
                else
                    j = 0;
                end
                
                while j<=numB&&elemB~=elemA
                    
                    j = j + 1;
                    elemB = epatchB(j);
                    
                end
                
                if elemA==elemB||i==numA
                    notdone = 0;
                end
                
            end
            
            if elemA==elemB % element interface
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
            nodeC = NodesOnElementCG(elemA,i);
            while i<5&&nodeC~=nodeA
                i = i + 1;
                nodeC = NodesOnElementCG(elemA,i);
            end
            nodeC = i;
            i = 1;
            nodeD = NodesOnElementCG(elemA,i);
            while i<5&&nodeD~=nodeB
                i = i + 1;
                nodeD = NodesOnElementCG(elemA,i);
            end
            nodeD = i;
            
            if nodeC*nodeD==2
                edgeA = 1;
                if nodeC>nodeD
                    nodeA = nodeC;
                    nodeC = nodeD;
                    nodeD = nodeA;
                end
            elseif nodeC*nodeD==6
                edgeA = 2;
                if nodeC>nodeD
                    nodeA = nodeC;
                    nodeC = nodeD;
                    nodeD = nodeA;
                end
            elseif nodeC*nodeD==12
                edgeA = 3;
                if nodeC>nodeD
                    nodeA = nodeC;
                    nodeC = nodeD;
                    nodeD = nodeA;
                end
            elseif nodeC*nodeD==3
                edgeA = 3;
                if nodeC<nodeD
                    nodeA = nodeC;
                    nodeC = nodeD;
                    nodeD = nodeA;
                end
            else % nodeC*nodeD==4
                edgeA = 4;
                if nodeC<nodeD
                    nodeA = nodeC;
                    nodeC = nodeD;
                    nodeD = nodeA;
                end
            end
            
            edgeB = elemedge(elemA,edgeA);
            
            if edgeB == 0
                
                numCL = numCL + 1;
                elemedge(elem,edge) = numCL;
                elemedge(elemA,edgeA) = numCL;
                SurfacesI(numCL,1) = NodesOnElementDG(elem,nloop(edge,2));
                SurfacesI(numCL,2) = NodesOnElementDG(elem,nloop(edge,1));
                SurfacesI(numCL,3) = NodesOnElementDG(elemA,nloop(edgeA,1));
                SurfacesI(numCL,4) = NodesOnElementDG(elemA,nloop(edgeA,2));
                SurfacesI(numCL,5) = elemA;
                SurfacesI(numCL,6) = elem;
                SurfacesI(numCL,7) = edgeA;
                SurfacesI(numCL,8) = edge;
                
            end
            
        else %domain boundary, add to SurfaceF
            
            numFL = numFL + 1;
            elemedge(elem,edge) = numFL;
            SurfacesF(numFL,1) = NodesOnElementDG(elem,nloop(edge,2));
            SurfacesF(numFL,2) = NodesOnElementDG(elem,nloop(edge,1));
            SurfacesF(numFL,3) = elem;
            SurfacesF(numFL,4) = edge;
            
        end
    
    end
    
end

clear('numD') % remove so that NL_FEA_Program doesn't find it
SurfacesI = SurfacesI(1:numCL,:);
SurfacesF = SurfacesF(1:numFL,:);

% Update boundary conditions
if exist('numBC','var') && numBC > 0
    NodeBCCG = NodeBC;
    numBCCG = numBC;
    [NodeBC,numBC] = updateCGDGnodelist(maxel,epnum,NodeCGDG,NodeBCCG,numBCCG);
end
if exist('numBCnp','var') && numBCnp > 0
    NodeBCnpCG = NodeBCnp;
    numBCnpCG = numBCnp;
    [NodeBCnp,numBCnp] = updateCGDGnodelist(maxel,epnum,NodeCGDG,NodeBCnpCG,numBCnpCG);
end
% Update loads
if exist('numNodalF','var') && numNodalF > 0
    NodeLoadCG = NodeLoad;
    numNodalFCG = numNodalF;
    [NodeLoad,numNodalF] = updateCGDGnodelist(maxel,epnum,NodeCGDG,NodeLoadCG,numNodalFCG);
end
if exist('numNodalFnp','var') && numNodalFnp > 0
    NodeLoadnpCG = NodeLoadnp;
    numNodalFnpCG = numNodalFnp;
    [NodeLoadnp,numNodalFnp] = updateCGDGnodelist(maxel,epnum,NodeCGDG,NodeLoadnpCG,numNodalFnpCG);
end

if exist('numSL','var') && numSL > 0
    SurfacesLCG = SurfacesL;
    numSLCG = numSL;
    [SurfacesL,numSL] = updateCGDGelemlist(ndm,nen,NodesOnElementCG,NodesOnElementDG,SurfacesLCG,numSLCG);
end
if exist('numSLnp','var') && numSLnp > 0
    SurfacesLnpCG = SurfacesLnp;
    numSLCGnp = numSLnp;
    [SurfacesLnp,numSLnp] = updateCGDGelemlist(ndm,nen,NodesOnElementCG,NodesOnElementDG,SurfacesLnpCG,numSLCGnp);
end
if exist('numBD','var') && numBD > 0
    SurfacesDCG = SurfacesD;
    numBDCG = numBD;
    [SurfacesD,numBD] = updateCGDGelemlist(ndm,nen,NodesOnElementCG,NodesOnElementDG,SurfacesDCG,numBDCG);
end
if exist('numUSL','var') && numUSL > 0
    SurfacesUCG = SurfacesU;
    numUSLCG = numUSL;
    [SurfacesU,numUSL] = updateCGDGelemlist(ndm,nen,NodesOnElementCG,NodesOnElementDG,SurfacesUCG,numUSLCG);
end
