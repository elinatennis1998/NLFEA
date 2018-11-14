function [SurfacesI,numSI] = getInterHoriz(Coordinates,nodeendR,nodeAR,nodeAL, ...
                             elemR,elemL,nodeincR,nodeincL,elemincR,elemincL,edgeR,edgeL)
%n1,n2,n3,n4,m1,m2,m3,m4,a,b,c,d,e,f)

% Function to specify nonconforming interface for two meshes with a
% horizontal interface; different element types are permitted

% nodeAR is the first node on the bottom surface
% nodeBR is the first node on the top surface (neither have to be exactly
%        at the start of the interface)
% elemR is the element for which nodeAR is in the upper left
% elemL is the element for which nodeBL is in the lower left
% nodeincR is the increment to nodes from one end of the element to the
%       next on the bottom surface
% elemincR is the increment to elements along the interface on the bottom

numSI = 1000;
SurfacesI = zeros(numSI,8);
numSI = 0;

Coordinates2 = Coordinates;
toler = 1e-12;

% Initialize
% nodeAR = (n1+1)*m1+1;
nodeBR = nodeAR + nodeincR;
nodeARx = Coordinates2(nodeAR,1);
nodeBRx = Coordinates2(nodeBR,1);
% nodeAL = b+1;
nodeBL = nodeAL+nodeincL;
nodeALx = Coordinates2(nodeAL,1);
nodeBLx = Coordinates2(nodeBL,1);
% elemL = a+1;
% elemR = n1*(m1-1)+1;

if nodeALx - nodeBRx > -toler  % nodeBRx <= nodeALx
while nodeALx - nodeBRx > -toler % nodeBRx < nodeALx
    nodeAR = nodeAR + nodeincR;
    nodeBR = nodeBR + nodeincR;
    elemR = elemR + elemincR;
    nodeBRx = Coordinates2(nodeBR,1);
end
elseif nodeARx - nodeBLx > -toler % nodeBLx <= nodeARx
while nodeARx - nodeBLx > -toler % nodeBLx < nodeARx
    nodeAL = nodeAL + nodeincL;
    nodeBL = nodeBL + nodeincL;
    elemL = elemL + elemincL;
    nodeBLx = Coordinates2(nodeBL,1);
end
end

nodeARx = Coordinates2(nodeAR,1);
nodeALx = Coordinates2(nodeAL,1);
nodeBRx = Coordinates2(nodeBR,1);
nodeBLx = Coordinates2(nodeBL,1);

while nodeBR <= nodeendR
    
    numSI = numSI + 1;
    if nodeALx - nodeARx > toler % nodeARx < nodeALx
        nodear = -1;
        nodeal = nodeAL;
    elseif nodeARx - nodeALx > toler % nodeARx > nodeALx
        nodear = nodeAR;
        nodeal = -1;
    else
        nodear = nodeAR;
        nodeal = nodeAL;
    end
    if nodeBRx - nodeBLx > toler % nodeBRx > nodeBLx
        nodebr = -1;
        nodebl = nodeBL;
    elseif nodeBLx - nodeBRx > toler % nodeBRx < nodeBLx
        nodebr = nodeBR;
        nodebl = -1;
    else
        nodebr = nodeBR;
        nodebl = nodeBL;
    end
    SurfacesI(numSI,:) = [nodear nodebr nodeal nodebl elemL elemR edgeL edgeR];
    
    if nodeBLx - nodeBRx > toler % nodeBRx < nodeBLx
        nodeAR = nodeBR;
        nodeBR = nodeBR + nodeincR;
        elemR = elemR + elemincR;
        nodeARx = nodeBRx;
        nodeBRx = Coordinates2(nodeBR,1);
    elseif nodeBRx - nodeBLx > toler % nodeBRx > nodeBLx
        nodeAL = nodeBL;
        nodeBL = nodeBL + nodeincL;
        elemL = elemL + elemincL;
        nodeALx = nodeBLx;
        nodeBLx = Coordinates2(nodeBL,1);
    else
        nodeAR = nodeBR;
        nodeBR = nodeBR + nodeincR;
        elemR = elemR + elemincR;
        nodeARx = nodeBRx;
        nodeBRx = Coordinates2(nodeBR,1);
        nodeAL = nodeBL;
        nodeBL = nodeBL + nodeincL;
        elemL = elemL + elemincL;
        nodeALx = nodeBLx;
        nodeBLx = Coordinates2(nodeBL,1);
    end
    
end


   
SurfacesI = SurfacesI(1:numSI,:);
