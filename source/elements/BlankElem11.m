function ElemE = BlankElem11(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,numEn)
        
% Purpose: Compute error norms for element
% Called by: Explicit_Error_Estimation.m
% Notes: General routine for assembling scalar quantities across all
%        elements, typically the element error norms but other quantities
%        could be used instead. The usual ordering is by derivative and
%        then by dof; see L_Elem3_2dVMS.m for an example.
% Example: L_Elem3_2dVMS.m

        ElemE = zeros(numEn,1);