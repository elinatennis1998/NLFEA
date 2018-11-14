function ElemS = BlankElem26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen,nestr)
        
% Purpose: Compute element stresses, e.g. at integration points
% Called by: FormFE.m
% Notes: Number of stresses set as nestr in pmatin; should be e.g. number
%        of integration points times number of stresses per point
% Example: NL_Elem2_2dM.m