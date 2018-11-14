function [ElemK,hr] = BlankElem21(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute stiffness matrix for element
% Called by: FormFE.m
% Notes: The ordering of dofs in the stiffness are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
%        For DG elements, the stiffness of the left and right sides may be
%        computed separately and then combined into the total stiffness
%        ElemK before exiting.
% Example: NL_Elem5_2dNSCST.m
        
        ElemK = zeros(nst);