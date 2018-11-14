function [ElemF,hr] = BlankElem06(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute residual vector for element
% Called by: FormFE.m
% Notes: The ordering of dofs in the force are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        The internal force vector should be the negative of its proper
%        value, e.g. -B'*sigma, due to the definition that the
%        Residual = F_external - F_internal.
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
%        For DG elements, the force of the left and right sides may be
%        computed separately and then combined into the total stiffness
%        ElemF before exiting.
% Example: NL_Elem5_2dNSCST.m
        
        ElemF = zeros(ndf*numel,1);