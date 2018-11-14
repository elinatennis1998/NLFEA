% c$Id:$
%       subroutine umacr6(lct,ctl)
% 
% c      * * F E A P * * A Finite Element Analysis Program
% 
% c....  Copyright (c) 1984-2014: Regents of the University of California
% c                               All rights reserved
% 
% c-----[--.----+----.----+----.-----------------------------------------]
% c     Modification log                                Date (dd/mm/year)
% c     Author: R.L. Taylor                                      8/14/2006
% c       Original version                                    01/11/2006
% c-----[--.----+----.----+----.-----------------------------------------]
% c      Purpose:  Convert each 8-node quadrilateral into a 10-node quad.
% 
% c      Inputs:
% c         lct       - Command character parameters
% c         ctl(3)    - Command numerical parameters
% 
% c      Outputs:
% c         File:     - 'Coord_10' : file containing nodal coordinates
% c                     for added node 10 values.
% c                   - 'Elemt_10' : file containing nodal connections
% c                     for added node 10 of elements that originally had
% c                     8 nodes.
% c-----[--.----+----.----+----.-----------------------------------------]
function [NodesOnElement,RegionOnElement,x,numnp,nen] = Q8toQ9(NodesOnElement,RegionOnElement,x,numel,ndm,nen)

    if nen ~= 8
        error('nen must be 8')
    end
    numnp = size(x,2);
    x = [x(1:ndm,1:numnp) zeros(ndm,numel)];
    NodesOnElement = [NodesOnElement(1:nen,1:numel); zeros(1,numel)];
    nen = 9;
    
    % Loop to add new node IDs into mesh and to assign them into NodesOnElement array
    % and to compute there averaged coordinates
    for elem = 1:numel
        nel = nnz(NodesOnElement(1:nen,elem));
        if nel == 8 % add new node to Q8 element
            numnp = numnp + 1;
            NodesOnElement(nen,elem) = numnp;
            xl = x(1:ndm,NodesOnElement(1:nen,elem));
            xmid = 0.50d0*(xl(:,5) + xl(:,6) + xl(:,7) + xl(:,8)) ...
                    - 0.25d0*(xl(:,1) + xl(:,2) + xl(:,3) + xl(:,4));
            x(1:ndm,numnp) = xmid;
        end
    end

end
