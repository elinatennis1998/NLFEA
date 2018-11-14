% Assemble Element Stiffness and Force into Model Stiffness and Force
%
% Tim Truster
% 10/2009
% UIUC

% Assembles assuming constraints are present
if numComp > 0
    
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
    Fd2(EGDOFTac) = Fd2(EGDOFTac) + ElemFTemp;
    
    ElemKTemp1 = zeros(lEGDOFTa ,lEGDOFTac);
    ElemKTemp2 = zeros(lEGDOFTac,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemKTemp1(:,EGDOFTan(i)) = ElemKTemp1(:,EGDOFTan(i)) + ElemK(ELDOFTa,ELDOFTa(i));
    end
    for i = 1:lEGDOFTa
        ElemKTemp2(EGDOFTan(i),:) = ElemKTemp2(EGDOFTan(i),:) + ElemKTemp1(i,:);
    end
    Kdd22(EGDOFTac,EGDOFTac) = Kdd22(EGDOFTac,EGDOFTac) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTa ,lEGDOFTic);
    ElemKTemp2 = zeros(lEGDOFTac,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemKTemp1(:,EGDOFTin(i)) = ElemKTemp1(:,EGDOFTin(i)) + ElemK(ELDOFTa,ELDOFTi(i));
    end
    for i = 1:lEGDOFTa
        ElemKTemp2(EGDOFTan(i),:) = ElemKTemp2(EGDOFTan(i),:) + ElemKTemp1(i,:);
    end
    Kdf2(EGDOFTac,EGDOFTic) = Kdf2(EGDOFTac,EGDOFTic) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTi ,lEGDOFTac);
    ElemKTemp2 = zeros(lEGDOFTic,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemKTemp1(:,EGDOFTan(i)) = ElemKTemp1(:,EGDOFTan(i)) + ElemK(ELDOFTi,ELDOFTa(i));
    end
    for i = 1:lEGDOFTi
        ElemKTemp2(EGDOFTin(i),:) = ElemKTemp2(EGDOFTin(i),:) + ElemKTemp1(i,:);
    end
    Kfd2(EGDOFTic,EGDOFTac) = Kfd2(EGDOFTic,EGDOFTac) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTi ,lEGDOFTic);
    ElemKTemp2 = zeros(lEGDOFTic,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemKTemp1(:,EGDOFTin(i)) = ElemKTemp1(:,EGDOFTin(i)) + ElemK(ELDOFTi,ELDOFTi(i));
    end
    for i = 1:lEGDOFTi
        ElemKTemp2(EGDOFTin(i),:) = ElemKTemp2(EGDOFTin(i),:) + ElemKTemp1(i,:);
    end
    Kff2(EGDOFTic,EGDOFTic) = Kff2(EGDOFTic,EGDOFTic) + ElemKTemp2;
    
else
    
    Fd2(EGDOFTa) = Fd2(EGDOFTa) + ElemF(ELDOFTa);
    Kdd22(EGDOFTa,EGDOFTa) = Kdd22(EGDOFTa,EGDOFTa) + ElemK(ELDOFTa,ELDOFTa);
    Kdf2(EGDOFTa,EGDOFTi) = Kdf2(EGDOFTa,EGDOFTi) + ElemK(ELDOFTa,ELDOFTi);
    Kfd2(EGDOFTi,EGDOFTa) = Kfd2(EGDOFTi,EGDOFTa) + ElemK(ELDOFTi,ELDOFTa);
    Kff2(EGDOFTi,EGDOFTi) = Kff2(EGDOFTi,EGDOFTi) + ElemK(ELDOFTi,ELDOFTi);

end
