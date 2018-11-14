% Tim Truster
% 06/15/2011
%
% Script to convert CG mesh (e.g. from block program output) into DG mesh
%
% Modified 07/23-24/2013 to not separate elements but to instead make
% interface elements for pressure discontinuity only

maxel = 12;
ndm = 2;
epatch = zeros(maxel,numnp);
epnum = zeros(numnp,1);
ixCG = ix;
numz = 0;

% Form epatch, epnum, DG nodes and connectivity

for elem = 1:numel

    nel = nnz(ixCG(elem,1:nen));
    
    for k = 1:nel % Loop over local Nodes

        node = ixCG(elem,k);

%   Add element to star list, increment number of elem in star

        q = epnum(node) + 1;
        epatch(q,node) = elem;		%epatch(nel,numnp)
        epnum(node) = q;			%epnum(numnp)


    end %k
    
end %j

numCL = 0;
numComp = 0;
SurfacesI = zeros(numel*4,8);
nloop = zeros(4,2);
elemedge = zeros(numel,4);

% Form SurfaceI, SurfaceL

for elem = 1:numel
    
    nel = nnz(ixCG(elem,1:nen));
    
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

        nodeA = ixCG(elem,nloop(edge,1));
        nodeB = ixCG(elem,nloop(edge,2));
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
            nodeC = ixCG(elemA,i);
            while i<5&&nodeC~=nodeA
                i = i + 1;
                nodeC = ixCG(elemA,i);
            end
            nodeC = i;
            i = 1;
            nodeD = ixCG(elemA,i);
            while i<5&&nodeD~=nodeB
                i = i + 1;
                nodeD = ixCG(elemA,i);
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
                SurfacesI(numCL,1) = ixCG(elem,nloop(edge,2));
                SurfacesI(numCL,2) = ixCG(elem,nloop(edge,1));
                SurfacesI(numCL,3) = ixCG(elemA,nloop(edgeA,1));
                SurfacesI(numCL,4) = ixCG(elemA,nloop(edgeA,2));
                SurfacesI(numCL,5) = elemA;
                SurfacesI(numCL,6) = elem;
                SurfacesI(numCL,7) = edgeA;
                SurfacesI(numCL,8) = edge;
                
            end
            
        end
    
    end
    
end

SurfacesI = SurfacesI(1:numCL,:);
