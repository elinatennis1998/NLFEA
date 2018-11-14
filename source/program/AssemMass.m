% Assemble Element Mass into Model Mass
%
% Tim Truster
% 8/2009
% UIUC

if lumping == 1
    ElemM = diag(sum(ElemM));
end

% Assembles assuming constraints are present
% if numComp > 0
    
    [EGDOFTac, EGDOFTam, EGDOFTan] = unique(EGDOFTa,'first');
    [EGDOFTic, EGDOFTim, EGDOFTin] = unique(EGDOFTi,'first');
    lEGDOFTac = length(EGDOFTac);
    lEGDOFTic = length(EGDOFTic);
    lEGDOFTa  = length(EGDOFTa );
    lEGDOFTi  = length(EGDOFTi );
    
    ElemMTemp1 = zeros(lEGDOFTa ,lEGDOFTac);
    ElemMTemp2 = zeros(lEGDOFTac,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemMTemp1(:,EGDOFTan(i)) = ElemMTemp1(:,EGDOFTan(i)) + ElemM(ELDOFTa,ELDOFTa(i));
    end
    for i = 1:lEGDOFTa
        ElemMTemp2(EGDOFTan(i),:) = ElemMTemp2(EGDOFTan(i),:) + ElemMTemp1(i,:);
    end
    Mdd11(EGDOFTac,EGDOFTac) = Mdd11(EGDOFTac,EGDOFTac) + ElemMTemp2;
    
    ElemMTemp1 = zeros(lEGDOFTa ,lEGDOFTic);
    ElemMTemp2 = zeros(lEGDOFTac,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemMTemp1(:,EGDOFTin(i)) = ElemMTemp1(:,EGDOFTin(i)) + ElemM(ELDOFTa,ELDOFTi(i));
    end
    for i = 1:lEGDOFTa
        ElemMTemp2(EGDOFTan(i),:) = ElemMTemp2(EGDOFTan(i),:) + ElemMTemp1(i,:);
    end
    Mdf1(EGDOFTac,EGDOFTic) = Mdf1(EGDOFTac,EGDOFTic) + ElemMTemp2;
    
    ElemMTemp1 = zeros(lEGDOFTi ,lEGDOFTac);
    ElemMTemp2 = zeros(lEGDOFTic,lEGDOFTac);
    for i = 1:lEGDOFTa
        ElemMTemp1(:,EGDOFTan(i)) = ElemMTemp1(:,EGDOFTan(i)) + ElemM(ELDOFTi,ELDOFTa(i));
    end
    for i = 1:lEGDOFTi
        ElemMTemp2(EGDOFTin(i),:) = ElemMTemp2(EGDOFTin(i),:) + ElemMTemp1(i,:);
    end
    Mfd(EGDOFTic,EGDOFTac) = Mfd(EGDOFTic,EGDOFTac) + ElemMTemp2;
    
    ElemMTemp1 = zeros(lEGDOFTi ,lEGDOFTic);
    ElemMTemp2 = zeros(lEGDOFTic,lEGDOFTic);
    for i = 1:lEGDOFTi
        ElemMTemp1(:,EGDOFTin(i)) = ElemMTemp1(:,EGDOFTin(i)) + ElemM(ELDOFTi,ELDOFTi(i));
    end
    for i = 1:lEGDOFTi
        ElemMTemp2(EGDOFTin(i),:) = ElemMTemp2(EGDOFTin(i),:) + ElemMTemp1(i,:);
    end
    Mff(EGDOFTic,EGDOFTic) = Mff(EGDOFTic,EGDOFTic) + ElemMTemp2;
    
% else
%     
%     Mdd11(EGDOFTa,EGDOFTa) = Mdd11(EGDOFTa,EGDOFTa) + ElemM(ELDOFTa,ELDOFTa);
%     Mdf1(EGDOFTa,EGDOFTi) = Mdf1(EGDOFTa,EGDOFTi) + ElemM(ELDOFTa,ELDOFTi);
%     Mfd(EGDOFTi,EGDOFTa) = Mfd(EGDOFTi,EGDOFTa) + ElemM(ELDOFTi,ELDOFTa);
%     Mff(EGDOFTi,EGDOFTi) = Mff(EGDOFTi,EGDOFTi) + ElemM(ELDOFTi,ELDOFTi);
% 
% end
