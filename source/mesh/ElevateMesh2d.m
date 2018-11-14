
maxel = 12;
epatch = zeros(maxel,numnp); % Array of elements connected to nodes
epnum = zeros(numnp,1);
Coordinates3 = zeros(numel*nen,ndm);
Coordinates3(1:numnp,:) = Coordinates;
NodesOnElementDG = zeros(numel,nen);
RegionOnElementDG = zeros(numel,1);
NodesOnElementCG = NodesOnElement;
RegionOnElementCG = RegionOnElement;
numz = 0;

% Form epatch, epnum, DG nodes and connectivity

for elem = 1:numel

    nel = nnz(ixCG(elem,1:nen));
    RegionOnElementDG(nen1) = RegionOnElement(nen1);
    
    for k = 1:nel % Loop over local Nodes

        node = NodesOnElementCG(elem,k);

%   Add element to star list, increment number of elem in star

        q = epnum(node) + 1;
        epatch(q,node) = elem;		%epatch(nel,numnp)
        epnum(node) = q;			%epnum(numnp)


    end %k
    
end %j

NodesOnElementDG = [NodesOnElementCG(:,1:nen) zeros(numel,nen)];
RegionOnElementDG = RegionOnElementCG;

numCL = 0;
numFL = 0;
numComp = 0;
SurfacesI = zeros(numel*4,8);
SurfacesF = zeros(numel,7); % All exposed faces
nloop = zeros(4,2);
elemedge = zeros(numel,4);

new_node = numnp;

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
                
                new_node = new_node + 1; % add a new node 
                nodeA = NodesOnElementDG(elem,nloop(edge,1));
                nodeB = NodesOnElementDG(elem,nloop(edge,2));
                nodeXY = (Coordinates(nodeA,:)+Coordinates(nodeB,:))/2;
                Coordinates3(new_node,:) = nodeXY;
                NodesOnElementDG(elem,nel+edge) = new_node; % put it in the slot for the midside node that it is, whether T or Q
                NodesOnElementDG(elemA,nel+edgeA) = new_node; % put it in the slot for the midside node that it is, whether T or Q
                numCL = numCL + 1;
                elemedge(elem,edge) = numCL;
                elemedge(elemA,edgeA) = numCL;
                
            end
            
        else %domain boundary, add to SurfaceF
            
            new_node = new_node + 1; % add a new node 
            nodeA = NodesOnElementDG(elem,nloop(edge,1));
            nodeB = NodesOnElementDG(elem,nloop(edge,2));
            nodeXY = (Coordinates(nodeA,:)+Coordinates(nodeB,:))/2;
            Coordinates3(new_node,:) = nodeXY;
            NodesOnElementDG(elem,nel+edge) = new_node; % put it in the slot for the midside node that it is, whether T or Q
            numFL = numFL + 1;
            elemedge(elem,edge) = numFL;
            
        end
    
    end
    
end

numnp = new_node;
Coordinates = Coordinates3(1:numnp,:);
NodesOnElement = NodesOnElementDG;
RegionOnElement = RegionOnElementDG;
if nen == 3
    nen = 6;
else
    nen = 9;
end
