% Assemble Element Stiffness and Force into Linear Model Stiffness and Force
%
% Tim Truster
% 04/17/2013
% UIUC

% Assembles assuming constraints are present
% if numComp > 0
    
    [EGDOFTac, EGDOFTam, EGDOFTan] = unique(EGDOFTa,'first');
    [EGDOFTic, EGDOFTim, EGDOFTin] = unique(EGDOFTi,'first');
    lEGDOFTac = length(EGDOFTac);
    lEGDOFTic = length(EGDOFTic);
    lEGDOFTa  = length(EGDOFTa );
    lEGDOFTi  = length(EGDOFTi );
    
    ElemKTemp1 = zeros(lEGDOFTa ,lEGDOFTac);
    ElemKTemp2 = zeros(lEGDOFTac,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemKTemp1(:,EGDOFTan(i)) = ElemKTemp1(:,EGDOFTan(i)) + ElemK(ELDOFTa,ELDOFTa(i));
    end
    for i = 1:lEGDOFTa
        ElemKTemp2(EGDOFTan(i),:) = ElemKTemp2(EGDOFTan(i),:) + ElemKTemp1(i,:);
    end
    KddLL(EGDOFTac,EGDOFTac) = KddLL(EGDOFTac,EGDOFTac) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTa ,lEGDOFTic);
    ElemKTemp2 = zeros(lEGDOFTac,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemKTemp1(:,EGDOFTin(i)) = ElemKTemp1(:,EGDOFTin(i)) + ElemK(ELDOFTa,ELDOFTi(i));
    end
    for i = 1:lEGDOFTa
        ElemKTemp2(EGDOFTan(i),:) = ElemKTemp2(EGDOFTan(i),:) + ElemKTemp1(i,:);
    end
    KdfLL(EGDOFTac,EGDOFTic) = KdfLL(EGDOFTac,EGDOFTic) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTi ,lEGDOFTac);
    ElemKTemp2 = zeros(lEGDOFTic,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemKTemp1(:,EGDOFTan(i)) = ElemKTemp1(:,EGDOFTan(i)) + ElemK(ELDOFTi,ELDOFTa(i));
    end
    for i = 1:lEGDOFTi
        ElemKTemp2(EGDOFTin(i),:) = ElemKTemp2(EGDOFTin(i),:) + ElemKTemp1(i,:);
    end
    KfdLL(EGDOFTic,EGDOFTac) = KfdLL(EGDOFTic,EGDOFTac) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTi ,lEGDOFTic);
    ElemKTemp2 = zeros(lEGDOFTic,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemKTemp1(:,EGDOFTin(i)) = ElemKTemp1(:,EGDOFTin(i)) + ElemK(ELDOFTi,ELDOFTi(i));
    end
    for i = 1:lEGDOFTi
        ElemKTemp2(EGDOFTin(i),:) = ElemKTemp2(EGDOFTin(i),:) + ElemKTemp1(i,:);
    end
    KffLL(EGDOFTic,EGDOFTic) = KffLL(EGDOFTic,EGDOFTic) + ElemKTemp2;
    
% else
%     
%     KddLL(EGDOFTa,EGDOFTa) = KddLL(EGDOFTa,EGDOFTa) + ElemK(ELDOFTa,ELDOFTa);
%     KdfLL(EGDOFTa,EGDOFTi) = KdfLL(EGDOFTa,EGDOFTi) + ElemK(ELDOFTa,ELDOFTi);
%     KfdLL(EGDOFTi,EGDOFTa) = KfdLL(EGDOFTi,EGDOFTa) + ElemK(ELDOFTi,ELDOFTa);
%     KffLL(EGDOFTi,EGDOFTi) = KffLL(EGDOFTi,EGDOFTi) + ElemK(ELDOFTi,ELDOFTi);
% 
% end
