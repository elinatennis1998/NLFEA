% Data Structure Allocation for Interface Solution Fields
% isw=61 must be set outside this routine
%
% Tim Truster
% 03/25/2014

% numIQ should be set in the element subroutine, case 61

% Check if integration triangular segments are used
if exist('xintt','var') && exist('ixt','var')
    maxnumseg = max(SurfacesI(:,10)-SurfacesI(:,9));
else
    maxnumseg = 1;
end

if ndm == 2
    nenseg = 2;
elseif ndm == 3
    nenseg = 4;
elseif ndm == 1
    nenseg = 1;
end

Iix = zeros(numSI*maxnumseg,nenseg+1); % interface segment connectivity
ICoordinates = zeros(numSI*maxnumseg*nenseg,ndm); % interface nodes
IElemSeg = zeros(2,numSI); % table with each DG element, the first and last segment associated with it

segment = 0; % segment counter
segnode = 0; % interface node ID counter

% Loop to set up nodal coordinates and connectivity of interface segments
for elem = numel-numSI+1:numel
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    
    %Determine element size parameters
    nel = nnz(NodesOnElement(elem,1:nen));
    nelP = getnelP(nel,ndm,nelP3,nelP4,nelP6,nelP9);
    nst = nel*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';
    
    ElemRout
    
   end
    
  end % ma
    
end % elem

numnpI = segnode;
numelI = segment;

% Convert to triangle segments if in 3-D
if nenseg == 4
    if max(Iix(1:numelI,nenseg)) == 0
        nenseg = 3;
    end
end

% Shrink arrays
Iix = Iix(1:numelI,[1:nenseg nenseg+1]); % interface segment connectivity
ICoordinates = ICoordinates(1:numnpI,1:ndm); % interface nodes
InterQuant = zeros(numIQ,numnpI,datastep); % interface nodes