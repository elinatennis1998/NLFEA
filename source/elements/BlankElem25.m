function ElemS = BlankElem25(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nelS,nen,npstr)
        
% Purpose: Projection/lumping of stresses from integration points to nodes
% Called by: FormS2.m
% Notes: Weighting can be accomplished either through elementarea or simply
%        the number of elements connected to a node.
%        Ordering of stresses for 2D is:
%        ElemS = [s_xx s_yy s_xy von_mises s_1 s_3 hydrostatic area]
%        Ordering of stresses for 3D is:
%        ElemS = [s_xx s_yy s_zz s_xy s_yz s_xz von_mises s_1 s_2 s_3 hydrostatic area]
%        Number of stresses set as npstr in pmatin
% Example: NL_Elem5_2dNSCST.m
        
        ElemS = zeros(nel,npstr+1);