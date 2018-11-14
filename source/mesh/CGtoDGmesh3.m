% Tim Truster
% 06/15/2011
%
% Script to convert CG mesh (e.g. from block program output) into DG mesh

maxel = 8;
ndm = 3;
epatch = zeros(maxel,numnp); % Array of elements connected to nodes
epnum = zeros(numnp,1);
NodeCGDG = epatch; % Array of CG nodes and the corresponding DG nodes
Coordinates3 = zeros(numel*nen,ndm);
NodesOnElementDG = ones(numel,nen);
NodesOnElementCG = NodesOnElement;
RegionOnElementCG = RegionOnElement;
RegionOnElementDG = RegionOnElement;
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
        epatch(q,node) = elem;		%epatch(nel,numnp) = element number
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
nloop = zeros(8,2);
elemedge = zeros(numel,8);

% Form SurfaceI, SurfaceL

for elem = 1:numel
    
    nel = nnz(NodesOnElementCG(elem,1:nen));
    
    if (nel==4||nel==10)
        
        nume = 4;
        nloop=[2 3 4
               1 3 4
               1 2 4
               1 2 3];
        
    else
        nume = 6;
        nloop=[1 4 8 5
               2 3 7 6
               1 2 6 5
               4 3 7 8
               1 2 3 4
               5 6 7 8];
        
    end
    
    %     Loop over edges of element
    if (nel==4||nel==10)
        
    for edge = 1:nume

        nodeA = NodesOnElementCG(elem,nloop(edge,1));
        nodeB = NodesOnElementCG(elem,nloop(edge,2));
        nodeC = NodesOnElementCG(elem,nloop(edge,3));
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
            nodeAA = NodesOnElementCG(elemA,i);
            while i<10&&nodeAA~=nodeA
                i = i + 1;
                nodeAA = NodesOnElementCG(elemA,i);
            end
            nodeAA = i;
            i = 1;
            nodeBB = NodesOnElementCG(elemA,i);
            while i<10&&nodeBB~=nodeB
                i = i + 1;
                nodeBB = NodesOnElementCG(elemA,i);
            end
            nodeBB = i;
            i = 1;
            nodeCC = NodesOnElementCG(elemA,i);
            while i<10&&nodeCC~=nodeC
                i = i + 1;
                nodeCC = NodesOnElementCG(elemA,i);
            end
            nodeCC = i;
            i = 1;
            sABC=sort([nodeAA nodeBB nodeCC]);
           for ii = 1:size(nloop,1)
               if norm(sABC-nloop(ii,:))==0
                  edgeA = ii;
               end     
           end
            
            edgeB = elemedge(elemA,edgeA);
            
            if edgeB == 0
                
                numCL = numCL + 1;
                elemedge(elem,edge) = numCL;
                elemedge(elemA,edgeA) = numCL;
                SurfacesI(numCL,1) = 0;
                SurfacesI(numCL,2) = 0;
                SurfacesI(numCL,3) = 0;
                SurfacesI(numCL,4) = 0;

                SurfacesI(numCL,5) = elem;
                SurfacesI(numCL,6) = elemA;
                SurfacesI(numCL,7) = edge;
                SurfacesI(numCL,8) = edgeA;
                
            end
            
        else %domain boundary, add to SurfacesF
            
            numFL = numFL + 1;
            elemedge(elem,edge) = numFL;
            SurfacesF(numFL,1) = NodesOnElementDG(elem,nloop(edge,2));
            SurfacesF(numFL,2) = NodesOnElementDG(elem,nloop(edge,1));
            SurfacesF(numFL,3) = elem;
            SurfacesF(numFL,4) = edge;
            
        end
    
    end
        
    else
        
    for edge = 1:nume

        nodeA = NodesOnElementCG(elem,nloop(edge,1));
        nodeB = NodesOnElementCG(elem,nloop(edge,2));
        nodeC = NodesOnElementCG(elem,nloop(edge,3));
        nodeD = NodesOnElementCG(elem,nloop(edge,4));
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
            nodeAA = NodesOnElementCG(elemA,i);
            while i<9&&nodeAA~=nodeA
                i = i + 1;
                nodeAA = NodesOnElementCG(elemA,i);
            end
            nodeAA = i;
            i = 1;
            nodeBB = NodesOnElementCG(elemA,i);
            while i<9&&nodeBB~=nodeB
                i = i + 1;
                nodeBB = NodesOnElementCG(elemA,i);
            end
            nodeBB = i;
            i = 1;
            nodeCC = NodesOnElementCG(elemA,i);
            while i<9&&nodeCC~=nodeC
                i = i + 1;
                nodeCC = NodesOnElementCG(elemA,i);
            end
            nodeCC = i;
            i = 1;
            nodeDD = NodesOnElementCG(elemA,i);
            while i<9&&nodeDD~=nodeD
                i = i + 1;
                nodeDD = NodesOnElementCG(elemA,i);
            end
            nodeDD = i;
%             if nodeC*nodeD==2
%                 edgeA = 1;
%                 if nodeC>nodeD
%                     nodeA = nodeC;
%                     nodeC = nodeD;
%                     nodeD = nodeA;
%                 end
%             elseif nodeC*nodeD==6
%                 edgeA = 2;
%                 if nodeC>nodeD
%                     nodeA = nodeC;
%                     nodeC = nodeD;
%                     nodeD = nodeA;
%                 end
%             elseif nodeC*nodeD==12
%                 edgeA = 3;
%                 if nodeC>nodeD
%                     nodeA = nodeC;
%                     nodeC = nodeD;
%                     nodeD = nodeA;
%                 end
%             elseif nodeC*nodeD==3
%                 edgeA = 3;
%                 if nodeC<nodeD
%                     nodeA = nodeC;
%                     nodeC = nodeD;
%                     nodeD = nodeA;
%                 end
%             else % nodeC*nodeD==4
%                 edgeA = 4;
%                 if nodeC<nodeD
%                     nodeA = nodeC;
%                     nodeC = nodeD;
%                     nodeD = nodeA;
%                 end
%            
%             end
           for ii = 1:size(nloop,1)
               if norm([nodeAA nodeBB nodeCC nodeDD]-nloop(ii,:))==0
                  edgeA = ii;
               end     
           end
            
            edgeB = elemedge(elemA,edgeA);
            
            if edgeB == 0
                
                numCL = numCL + 1;
                elemedge(elem,edge) = numCL;
                elemedge(elemA,edgeA) = numCL;
%                 SurfacesI(numCL,1) = NodesOnElementDG(elem,nloop(edge,2));
%                 SurfacesI(numCL,2) = NodesOnElementDG(elem,nloop(edge,1));
%                 SurfacesI(numCL,3) = NodesOnElementDG(elemA,nloop(edgeA,1));
%                 SurfacesI(numCL,4) = NodesOnElementDG(elemA,nloop(edgeA,2));
                SurfacesI(numCL,1) = 0;
                SurfacesI(numCL,2) = 0;
                SurfacesI(numCL,3) = 0;
                SurfacesI(numCL,4) = 0;

                SurfacesI(numCL,5) = elem;
                SurfacesI(numCL,6) = elemA;
                SurfacesI(numCL,7) = edge;
                SurfacesI(numCL,8) = edgeA;
                
            end
            
        else %domain boundary, add to SurfacesF
            
            numFL = numFL + 1;
            elemedge(elem,edge) = numFL;
            SurfacesF(numFL,1) = NodesOnElementDG(elem,nloop(edge,2));
            SurfacesF(numFL,2) = NodesOnElementDG(elem,nloop(edge,1));
            SurfacesF(numFL,3) = elem;
            SurfacesF(numFL,4) = edge;
            
        end
    
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
