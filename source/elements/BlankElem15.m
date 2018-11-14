function ElemF = BlankElem15(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute body force for element
% Called by: pbodyf.m
% Notes: The ordering of dofs in the external load vector are as follows:
%        ElemF = [node1_dof1 node1_dof2 ... node1_dofndf node2_dof1 ...
%                 node2_dofndf node3_dof1 ... nodenen_dofndf]';
%        Body forces are treated as external loads, which can be scaled by
%        a proportional loading parameter.
%        ElemF is initialized as zero inside pbodyf.m and does not need to
%        be reinitialized here.
% Example: NL_Elem2_2dM.m