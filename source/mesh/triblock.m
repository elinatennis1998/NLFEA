function [x,NodesOnElement,RegionOnElement,numnp,numel] = triblock(type,ninc,node1,elmt1,mat,ntyp,xl,nen,x,NodesOnElement,RegionOnElement)

% Block routine to emulate behavior of FEAP meshing routine TRIBlock
% Handles 2D T3 elements

if strcmp(type,'cart') == 1

    if nargin < 8 || nargin == 9
        error('At least 8 arguments are required')
    elseif nargin == 8
        x = zeros(0,0);
        NodesOnElement = zeros(0,0);
    elseif nargin == 10
        %x and NodesOnElement are given
    end
    
    xll = zeros(2,6);
    ixl = zeros(6,1);
    for i = 1:length(xl)
        ind = xl(i,1);
        ixl(ind) = ind;
        xll(1,ind) = xl(i,2);
        xll(2,ind) = xl(i,3);
    end
    nr = ninc;

    nod1 = node1;
    nuel1 = elmt1;
    ma = mat;
    if ntyp == 0
        ntyp = 1;
    end
    ndm = 2;
%     nm = 6;
    nen1 = nen+1;

% %     Determine last element to be generated
% 
%       if(ntyp == 1) %then
%         nf = nuel1 + nr*nr - 1;
%       else
%         nf = nuel1 + nr*nr/4 - 1;
%       end
% 
% %     Determine last node number to be generated
% 
%       ng = nod1 + (nr+1)*(nr+2)/2 - 1;

%     Form block of elements

    [x,NodesOnElement,RegionOnElement] = trblk(nr,xll,ixl,x,NodesOnElement,RegionOnElement,ndm,nod1,nuel1,nen1,ma,ntyp);
    numnp = length(x);
    numel = length(NodesOnElement);
    
end