function ElemE = BlankElem12(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelP,nen)
        
% Purpose: Compute energy for element
% Called by: FormFE.m
% Notes: Typically used for energy-preserving dynamic algorithms. May
%        require only computing the potential energy or also the kinetic
%        energy; see the specific algorithm's implementation for
%        determining what is necessary.
% Example: NL_Elem3_2d.m

        ElemE = 0;