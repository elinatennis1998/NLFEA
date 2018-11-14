% Assemble Interface Stiffness and Force into Model Stiffness and Force
%
% Tim Truster
% 3/2011
% UIUC

if mnr == 0 || iter == 0 %assemble K only for full Newton or first step
    
% Assembles assuming constraints are present
if numComp > 0
    
    [EGDOFTaLc, EGDOFTaLm, EGDOFTaLn] = unique(EGDOFTaL,'first');
    [EGDOFTaRc, EGDOFTaRm, EGDOFTaRn] = unique(EGDOFTaR,'first');
    [EGDOFTiLc, EGDOFTiLm, EGDOFTiLn] = unique(EGDOFTiL,'first');
    [EGDOFTiRc, EGDOFTiRm, EGDOFTiRn] = unique(EGDOFTiR,'first');
    lEGDOFTaLc = length(EGDOFTaLc);
    lEGDOFTaRc = length(EGDOFTaRc);
    lEGDOFTaL  = length(EGDOFTaL);
    lEGDOFTaR  = length(EGDOFTaR);
    lEGDOFTiLc = length(EGDOFTiLc);
    lEGDOFTiRc = length(EGDOFTiRc);
    lEGDOFTiL  = length(EGDOFTiL);
    lEGDOFTiR  = length(EGDOFTiR);
    
    ElemKTemp1 = zeros(lEGDOFTaL ,lEGDOFTaLc);
    ElemKTemp2 = zeros(lEGDOFTaLc,lEGDOFTaLc);
    for i = 1:lEGDOFTaL
        ElemKTemp1(:,EGDOFTaLn(i)) = ElemKTemp1(:,EGDOFTaLn(i)) + ElemKLL(ELDOFTaL,ELDOFTaL(i));
    end
    for i = 1:lEGDOFTaL
        ElemKTemp2(EGDOFTaLn(i),:) = ElemKTemp2(EGDOFTaLn(i),:) + ElemKTemp1(i,:);
    end
    Kdd11(EGDOFTaLc,EGDOFTaLc) = Kdd11(EGDOFTaLc,EGDOFTaLc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaL ,lEGDOFTaRc);
    ElemKTemp2 = zeros(lEGDOFTaLc,lEGDOFTaRc);
    for i = 1:lEGDOFTaR
        ElemKTemp1(:,EGDOFTaRn(i)) = ElemKTemp1(:,EGDOFTaRn(i)) + ElemKLR(ELDOFTaL,ELDOFTaR(i));
    end
    for i = 1:lEGDOFTaL
        ElemKTemp2(EGDOFTaLn(i),:) = ElemKTemp2(EGDOFTaLn(i),:) + ElemKTemp1(i,:);
    end
    Kdd11(EGDOFTaLc,EGDOFTaRc) = Kdd11(EGDOFTaLc,EGDOFTaRc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaR ,lEGDOFTaLc);
    ElemKTemp2 = zeros(lEGDOFTaRc,lEGDOFTaLc);
    for i = 1:lEGDOFTaL
        ElemKTemp1(:,EGDOFTaLn(i)) = ElemKTemp1(:,EGDOFTaLn(i)) + ElemKRL(ELDOFTaR,ELDOFTaL(i));
    end
    for i = 1:lEGDOFTaR
        ElemKTemp2(EGDOFTaRn(i),:) = ElemKTemp2(EGDOFTaRn(i),:) + ElemKTemp1(i,:);
    end
    Kdd11(EGDOFTaRc,EGDOFTaLc) = Kdd11(EGDOFTaRc,EGDOFTaLc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaR ,lEGDOFTaRc);
    ElemKTemp2 = zeros(lEGDOFTaRc,lEGDOFTaRc);
    for i = 1:lEGDOFTaR
        ElemKTemp1(:,EGDOFTaRn(i)) = ElemKTemp1(:,EGDOFTaRn(i)) + ElemKRR(ELDOFTaR,ELDOFTaR(i));
    end
    for i = 1:lEGDOFTaR
        ElemKTemp2(EGDOFTaRn(i),:) = ElemKTemp2(EGDOFTaRn(i),:) + ElemKTemp1(i,:);
    end
    Kdd11(EGDOFTaRc,EGDOFTaRc) = Kdd11(EGDOFTaRc,EGDOFTaRc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaL ,lEGDOFTiLc);
    ElemKTemp2 = zeros(lEGDOFTaLc,lEGDOFTiLc);
    for i = 1:lEGDOFTiL
        ElemKTemp1(:,EGDOFTiLn(i)) = ElemKTemp1(:,EGDOFTiLn(i)) + ElemKLL(ELDOFTaL,ELDOFTiL(i));
    end
    for i = 1:lEGDOFTaL
        ElemKTemp2(EGDOFTaLn(i),:) = ElemKTemp2(EGDOFTaLn(i),:) + ElemKTemp1(i,:);
    end
    Kdf1(EGDOFTaLc,EGDOFTiLc) = Kdf1(EGDOFTaLc,EGDOFTiLc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaL ,lEGDOFTiRc);
    ElemKTemp2 = zeros(lEGDOFTaLc,lEGDOFTiRc);
    for i = 1:lEGDOFTiR
        ElemKTemp1(:,EGDOFTiRn(i)) = ElemKTemp1(:,EGDOFTiRn(i)) + ElemKLR(ELDOFTaL,ELDOFTiR(i));
    end
    for i = 1:lEGDOFTaL
        ElemKTemp2(EGDOFTaLn(i),:) = ElemKTemp2(EGDOFTaLn(i),:) + ElemKTemp1(i,:);
    end
    Kdf1(EGDOFTaLc,EGDOFTiRc) = Kdf1(EGDOFTaLc,EGDOFTiRc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaR ,lEGDOFTiLc);
    ElemKTemp2 = zeros(lEGDOFTaRc,lEGDOFTiLc);
    for i = 1:lEGDOFTiL
        ElemKTemp1(:,EGDOFTiLn(i)) = ElemKTemp1(:,EGDOFTiLn(i)) + ElemKRL(ELDOFTaR,ELDOFTiL(i));
    end
    for i = 1:lEGDOFTaR
        ElemKTemp2(EGDOFTaRn(i),:) = ElemKTemp2(EGDOFTaRn(i),:) + ElemKTemp1(i,:);
    end
    Kdf1(EGDOFTaRc,EGDOFTiLc) = Kdf1(EGDOFTaRc,EGDOFTiLc) + ElemKTemp2;
    
    ElemKTemp1 = zeros(lEGDOFTaR ,lEGDOFTiRc);
    ElemKTemp2 = zeros(lEGDOFTaRc,lEGDOFTiRc);
    for i = 1:lEGDOFTiR
        ElemKTemp1(:,EGDOFTiRn(i)) = ElemKTemp1(:,EGDOFTiRn(i)) + ElemKRR(ELDOFTaR,ELDOFTiR(i));
    end
    for i = 1:lEGDOFTaR
        ElemKTemp2(EGDOFTaRn(i),:) = ElemKTemp2(EGDOFTaRn(i),:) + ElemKTemp1(i,:);
    end
    Kdf1(EGDOFTaRc,EGDOFTiRc) = Kdf1(EGDOFTaRc,EGDOFTiRc) + ElemKTemp2;
    
else
    
    Kdd11(EGDOFTaL,EGDOFTaL) = Kdd11(EGDOFTaL,EGDOFTaL) + ElemKLL(ELDOFTaL,ELDOFTaL);
    Kdd11(EGDOFTaL,EGDOFTaR) = Kdd11(EGDOFTaL,EGDOFTaR) + ElemKLR(ELDOFTaL,ELDOFTaR);
    Kdd11(EGDOFTaR,EGDOFTaL) = Kdd11(EGDOFTaR,EGDOFTaL) + ElemKRL(ELDOFTaR,ELDOFTaL);
    Kdd11(EGDOFTaR,EGDOFTaR) = Kdd11(EGDOFTaR,EGDOFTaR) + ElemKRR(ELDOFTaR,ELDOFTaR);
    
    Kdf1(EGDOFTaL,EGDOFTiL) = Kdf1(EGDOFTaL,EGDOFTiL) + ElemKLL(ELDOFTaL,ELDOFTiL);
    Kdf1(EGDOFTaL,EGDOFTiR) = Kdf1(EGDOFTaL,EGDOFTiR) + ElemKLR(ELDOFTaL,ELDOFTiR);
    Kdf1(EGDOFTaR,EGDOFTiL) = Kdf1(EGDOFTaR,EGDOFTiL) + ElemKRL(ELDOFTaR,ELDOFTiL);
    Kdf1(EGDOFTaR,EGDOFTiR) = Kdf1(EGDOFTaR,EGDOFTiR) + ElemKRR(ELDOFTaR,ELDOFTiR);

end

% locind1 = 0;
% for ie1 = 1:nel2L
%     for l1 = 1:ndf
%         locind1 = locind1 + 1;
%         growL = EDOFTL(locind1);
%         if growL > 0
%             if growL <= neq1
%                 locind2 = 0;
%                 for ie2 = 1:nel2L
%                     for l2 = 1:ndf
%                         locind2 = locind2 + 1;
%                         gcolL = EDOFTL(locind2);
%                         if gcolL <= neq1 && gcolL > 0
%                             Kdd11(growL, gcolL) = Kdd11(growL, gcolL) + ElemKLL(locind1,locind2); %#ok<AGROW>
%                         end
%                     end %l2
%                 end %ie2
%                 locind2 = 0;
%                 for ie2 = 1:nel2R
%                     for l2 = 1:ndf
%                         locind2 = locind2 + 1;
%                         gcolR = EDOFTR(locind2);
%                         if gcolR <= neq1 && gcolR > 0
%                             Kdd11(growL, gcolR) = Kdd11(growL, gcolR) + ElemKLR(locind1,locind2); %#ok<AGROW>
%                         end
%                     end %l2
%                 end %ie2
%             end
%         end
%     end %l1
% end %ie1
% locind1 = 0;
% for ie1 = 1:nel2R
%     for l1 = 1:ndf
%         locind1 = locind1 + 1;
%         growR = EDOFTR(locind1);
%         if growR > 0
%             if growR <= neq1
%                 locind2 = 0;
%                 for ie2 = 1:nel2L
%                     for l2 = 1:ndf
%                         locind2 = locind2 + 1;
%                         gcolL = EDOFTL(locind2);
%                         if gcolL <= neq1 && gcolL > 0
%                             Kdd11(growR, gcolL) = Kdd11(growR, gcolL) + ElemKRL(locind1,locind2); %#ok<AGROW>
%                         end
%                     end %l2
%                 end %ie2
%                 locind2 = 0;
%                 for ie2 = 1:nel2R
%                     for l2 = 1:ndf
%                         locind2 = locind2 + 1;
%                         gcolR = EDOFTR(locind2);
%                         if gcolR <= neq1 && gcolR > 0
%                             Kdd11(growR, gcolR) = Kdd11(growR, gcolR) + ElemKRR(locind1,locind2); %#ok<AGROW>
%                         end
%                     end %l2
%                 end %ie2
%             end
%         end
%     end %l1
% end %ie1

end

if rflag == 1

    if numComp > 0
    ElemFTemp = zeros(lEGDOFTaLc,1);
    
    for i = 1:lEGDOFTaL
        ElemFTemp(EGDOFTaLn(i)) = ElemFTemp(EGDOFTaLn(i)) + ElemFL(ELDOFTaL(i));
    end
    Fd1(EGDOFTaLc) = Fd1(EGDOFTaLc) + ElemFTemp;

    ElemFTemp = zeros(lEGDOFTaRc,1);
    for i = 1:lEGDOFTaR
        ElemFTemp(EGDOFTaRn(i)) = ElemFTemp(EGDOFTaRn(i)) + ElemFR(ELDOFTaR(i));
    end
    Fd1(EGDOFTaRc) = Fd1(EGDOFTaRc) + ElemFTemp;
    
    else
        
    Fd1(EGDOFTaL) = Fd1(EGDOFTaL) + ElemFL(ELDOFTaL);
    Fd1(EGDOFTaR) = Fd1(EGDOFTaR) + ElemFR(ELDOFTaR);
    
    Fd3(EGDOFTiL) = Fd3(EGDOFTiL) + ElemFL(ELDOFTiL);
    Fd3(EGDOFTiR) = Fd3(EGDOFTiR) + ElemFR(ELDOFTiR);
    
    end

%     locind1 = 0;
%     for ie1 = 1:nel2L
%         for l1 = 1:ndf
%             locind1 = locind1 + 1;
%             growL = EDOFTL(locind1);
%             if growL > 0
%                 if growL <= neq1
%                     Fd1(growL) = Fd1(growL) + ElemFL(locind1); %#ok<AGROW>
%                 end
%             end
%         end %l1
%     end %ie1
%     locind1 = 0;
%     for ie1 = 1:nel2R
%         for l1 = 1:ndf
%             locind1 = locind1 + 1;
%             growR = EDOFTR(locind1);
%             if growR > 0
%                 if growR <= neq1
%                     Fd1(growR) = Fd1(growR) + ElemFR(locind1); %#ok<AGROW>
%                 end
%             end
%         end %l1
%     end %ie1

end



%     ElemFL = zeros(nstL,1);
%     ElemFR = zeros(nstR,1);
% 
%     %Get Gc
%     loc1 = 0;
%     for o = 1:nel2L
%         for l = 1:ndf
%             loc1 = loc1 + 1;
%             sum1 = 0;
%             loc2 = 0;
%             for i = 1:nel2L
%                 for k = 1:ndf
%                     loc2 = loc2 + 1;
%                     sum1 = sum1 + ElemKLL(loc1,loc2)*ulL(k,i);
%                 end
%             end
%             loc2 = 0;
%             for i = 1:nel2R
%                 for k = 1:ndf
%                     loc2 = loc2 + 1;
%                     sum1 = sum1 + ElemKLR(loc1,loc2)*ulR(k,i);
%                 end
%             end
%             ElemFL(loc1) = ElemFL(loc1) - sum1;
%         end
%     end
%     loc1 = 0;
%     for o = 1:nel2R
%         for l = 1:ndf
%             loc1 = loc1 + 1;
%             sum1 = 0;
%             loc2 = 0;
%             for i = 1:nel2L
%                 for k = 1:ndf
%                     loc2 = loc2 + 1;
%                     sum1 = sum1 + ElemKRL(loc1,loc2)*ulL(k,i);
%                 end
%             end
%             loc2 = 0;
%             for i = 1:nel2R
%                 for k = 1:ndf
%                     loc2 = loc2 + 1;
%                     sum1 = sum1 + ElemKRR(loc1,loc2)*ulR(k,i);
%                 end
%             end
%             ElemFR(loc1) = ElemFR(loc1) - sum1;
%         end
%     end