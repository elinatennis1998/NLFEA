% Assemble Element Force into Model Force, user surface loads
%
% Tim Truster
% 12/30/2013
% UIUC

% Assembles assuming constraints are present
if numComp > 0
    
    [EGDOFTac, EGDOFTam, EGDOFTan] = unique(EGDOFTa,'first');
    lEGDOFTac = length(EGDOFTac);
    lEGDOFTa  = length(EGDOFTa );

    ElemFTemp = zeros(lEGDOFTac,1);
    for i = 1:lEGDOFTa
        ElemFTemp(EGDOFTan(i)) = ElemFTemp(EGDOFTan(i)) + ElemF(ELDOFTa(i));
    end
    FcU(EGDOFTac) = FcU(EGDOFTac) + ElemFTemp;
    
else
    
    FcU(EGDOFTa) = FcU(EGDOFTa) + ElemF(ELDOFTa);

end