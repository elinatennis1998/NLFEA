% Assemble Element Mass into Linear Model Mass
%
% Tim Truster
% 2/11/2014
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
    MddLL(EGDOFTac,EGDOFTac) = MddLL(EGDOFTac,EGDOFTac) + ElemMTemp2;
    
% else
%     
%     MddLL(EGDOFTa,EGDOFTa) = MddLL(EGDOFTa,EGDOFTa) + ElemM(ELDOFTa,ELDOFTa);
% 
% end
