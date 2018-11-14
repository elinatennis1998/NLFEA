function [ElemM,hr] = BlankElem05(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute mass matrix for element
% Called by: FormFE.m
% Notes: The ordering of dofs in the mass are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
%        Mass matrices are only used by some dynamic algorithms, and then
%        possibly only to compute the initial accelarations. In other
%        cases, the mass is recomputed within a composite stiffness matrix
%        in the isw=3 task.
% Example: NL_Elem3_2d.m
        
        ElemM = zeros(nst);