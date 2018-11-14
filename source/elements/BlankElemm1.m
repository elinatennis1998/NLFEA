function ElemF = BlankElemm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute surface traction for an element
% Called by: ploadi.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Integration of the force can be handled in a standard fashion by
%        rotating the parameteric space of the current element to match a
%        template face over which integration is always performed.
% Example: L_Elem3_2dVMS.m
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Perform integration here
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);