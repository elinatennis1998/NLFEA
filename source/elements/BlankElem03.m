function [ElemK,ElemF,hr] = BlankElem03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute stiffness matrix and residual vector for element
% Called by: FormFE.m
% Notes: The ordering of dofs in the stiffness/force are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        The internal force vector should be the negative of its proper
%        value, e.g. -B'*sigma, due to the definition that the
%        Residual = F_external - F_internal.
%        While nst = ndf*nen in general, typically the assembly process
%        only requires the values for i=1:ndf*nel, where nel is the number
%        of nodes on the current element.
%        For DG elements, the stiffness of the left and right sides may be
%        computed separately and then combined into the total stiffness
%        ElemK before exiting.
% Example: L_Elem3_2dVMS.m
        
        ElemK = zeros(ndf*numel,ndf*numel);
        ElemF = zeros(ndf*numel,1);
        k1 = PatchE*PatchA/(L/3);
        k2 = PatchE*PatchA/(2*L/3);
        
        ElemeK = [k1 -k1 0
            -k1 k1+k2 -k2
            0 -k2 k2];
        P = 100;
        f1 = q*(L/3)/2.*[1;1]:
        f2 = q*(2*L/3)/2.*[1;1];
        ElemF = [f1;f1;0]+[0;f2;f2]+[0;0;P];
        
        