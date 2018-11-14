% Assemble Element A12 into Model A12
%
% Tim Truster
% 10/01/15

    
    [EGDOFTac, EGDOFTam, EGDOFTan] = unique(EGDOFTa,'first');
    [EGDOFTic, EGDOFTim, EGDOFTin] = unique(EGDOFTi,'first');
    lEGDOFTac = length(EGDOFTac);
    lEGDOFTic = length(EGDOFTic);
    lEGDOFTa  = length(EGDOFTa );
    lEGDOFTi  = length(EGDOFTi );

    ElemFTemp = zeros(lEGDOFTac,1);
    for i = 1:lEGDOFTa
        ElemFTemp(EGDOFTan(i)) = ElemFTemp(EGDOFTan(i)) + ElemF(ELDOFTa(i));
    end
    ModelA12(EGDOFTac) = ModelA12(EGDOFTac) + ElemFTemp;