% Assemble Quantities from Model Routine
%
% Tim Truster
% 7/2009
% UIUC
%
% NOTE 02/11/2014: for dynamic analyses with linear elements, the stiffness
% and mass matrix must be kept separate in order to shorten the calculation
% time. Thus, in the element routine mass is handled by isw=5 and stiffness
% by isw=3, and the% results are combined within the main program in Mdd11 
% and Kdd11.
% However, for nonlinear analyses the FEAP-standard can be used where
% K and M are combined together as the Mstar matrix. The code assumes that
% the entire static and dynamic contribution is computed by isw=3 and then
% assembled into Kdd11.
% This has been updated for transient = 1,2,5
% transient=3 still assumes separate matrices for M and K

% Initialize Model-level arrays
hflgu = 0;
h3flgu = 0;

switch isw%Task Switch
    
    case 1 %Get Material Properties
        
    case 3 %Get Stiffness, Force
        hflgu = 1;
        h3flgu = 1;
        if initializeLinKF && numberlinear > 0 %initialize linear stiffness
            KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
            KfdLL = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
            KffLL = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
            Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        elseif numberlinear > 0 %load linear stiffness
            Kdd11 = KddLL;
            Kdf1 = KdfLL;
        else %no linear stiffness, so recompute each time
            Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        end
        Kfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Kff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemStifForc';
        switch transient
            case {-1,0}
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Kdd11*ModelDx - Kdf1*gBC;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case {1,2}
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
                    Fd1 = Fd1 - Mdd11*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Kdd11 = coeffm*Mdd11 + coeffk*Kdd11;
                end
%                 F1n = F1n - alpha*Kdd11*ModelDx - alpha*Kdf1*gBC;
            case 3
                Fd1 = coeff1*Fext1;
            case 4
                Fd1 = coeff1*Fext1;
            case 5
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
                    Fd1 = Fd1 - Mdd11*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Kdd11 = coeffm*Mdd11 + coeffk*Kdd11;
                end
            case 6
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
                end
            case 8
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
                end
            case 10
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
                    Fd1 = Fd1 - Mdd11*ModelAx - Kdd11*ModelDx - Kdf1*gBC;
                    Kdd11 = Mdd11;
                end
            case 11
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Kdd11*ModelDx - Kdf1*gBC;
                end
                Fd3 = zeros(nieq,1);
            case 12
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
                    Fd1 = Fd1 - Mdd11*ModelAx - Kdd11*ModelDx - Kdf1*gBC;
                    Kdd11 = Mdd11;
                end
            case 13
                Fd1 = coeff1*Fext1;
            case 14
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
%                     Fd1 = Fd1 - Mdd11*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Fd1 = Fd1 - Mdd11*(coeffml2*ModelVx+coeffml1*ModelVxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Kdd11 = coeffm*Mdd11 + coeffk*Kdd11;
                end
        end
%             if numD > 0
%                 iswtemp = isw;
%                 isw = -4;
%                 FormDN
%                 isw = iswtemp;
%             end
    case 5 %Get Mass
        if initializeLinKF && numberlinear > 0 %initialize linear stiffness
            MddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            Mdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
        elseif numberlinear > 0 %load linear stiffness
            Mdd11 = MddLL;
        else %no linear stiffness, so recompute each time
            Mdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
        end
        Mdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        Mfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Mff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemMass';
    case 6 %Get Force
        hflgu = 1;
        h3flgu = 1;
        AssemQuant = 'AssemForc';
%         if ~initializeLinKF && numberlinear > 0 %load linear stiffness
%             Kdd11 = KddLL;
%             Kdf1 = KdfLL;
%         else %no linear stiffness, so recompute each time
%             Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
%             Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
%         end
        switch transient
            case {-1,0}
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - KddLL*ModelDx - KdfLL*gBC;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case {1,2}
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - MddLL*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - KddLL*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - KdfLL*(coeffkl2*gBC+coeffkl1*gBC_n);
                    KddLL = coeffm*MddLL + coeffk*KddLL;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
%                 F1n = F1n - alpha*KddLL*ModelDx - alpha*KdfLL*gBC;
            case 3
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
            case 4
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
            case 5
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - MddLL*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - KddLL*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - KdfLL*(coeffkl2*gBC+coeffkl1*gBC_n);
                    KddLL = coeffm*MddLL + coeffk*KddLL;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 6
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF % THIS WILL NOT WORK ANYMORE - NEED TO STORE K AND M in NLFEA somewhere at start
                    Fd1 = Fd1 - MddLL*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - KddLL*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - KdfLL*gBC;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 8
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF % THIS WILL NOT WORK ANYMORE - NEED TO STORE K AND M in NLFEA somewhere at start
                    Fd1 = Fd1 - MddLL*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - KddLL*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - KdfLL*gBC;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 10
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - MddLL*ModelAx - KddLL*ModelDx - KdfLL*gBC;
                    KddLL = MddLL;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 11
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - KddLL*ModelDx - KdfLL*gBC;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 12
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
                    Fd1 = Fd1 - Mdd11*ModelAx - Kdd11*ModelDx - Kdf1*gBC;
                    Kdd11 = Mdd11;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
            case 13
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
            case 14
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Mdd11 = MddLL;
%                     Fd1 = Fd1 - Mdd11*(coeffml2*ModelAx+coeffml1*ModelAxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Fd1 = Fd1 - Mdd11*(coeffml2*ModelVx+coeffml1*ModelVxn_1) - Kdd11*(coeffkl2*ModelDx+coeffkl1*ModelDxn_1) - Kdf1*(coeffkl2*gBC+coeffkl1*gBC_n);
                    Kdd11 = coeffm*Mdd11 + coeffk*Kdd11;
                    Fd3 = Fd3 - KfdLL*ModelDx - KffLL*gBC;
                end
        end
    case 9 %Get global error
        if transient > 0
            Fd1 = coeff1*Fext1;
            Fd3 = zeros(nieq,1);
            AssemQuant = 'AssemForc';
        else
            Fd1 = zeros(neq,1);
            Fd3 = zeros(nieq,1);
            AssemQuant = 'AssemForc';
        end
    case 11 %Get Error Indicators
        AssemQuant = 'AssemEner';
        Energy = zeros(numEn,1);
    case 12 % system energy
        SysEner = 0;
        AssemQuant = 'AssemEner2';
    case 13 % plastic dissipation
        PlasDiss = 0;
        AssemQuant = 'AssemDiss';
    case 21 %Get Stiffness
        hflgu = 1;
        h3flgu = 1;
        if initializeLinKF && numberlinear > 0 %initialize linear stiffness
            KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
            KfdLL = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
            KffLL = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
            Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        elseif numberlinear > 0 %load linear stiffness
            Kdd11 = KddLL;
            Kdf1 = KdfLL;
        else %no linear stiffness, so recompute each time
            Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        end
        Kfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Kff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemStif';
    case 22
        AssemQuant = 'AssemStre';
    case 26
        AssemQuant = 'AssemEStre';
    case 24
        AssemQuant = 'AssemPlas';
    case 40
        AssemQuant = '';
        hflgu = 1;
        h3flgu = 1;
    case 51
        AssemQuant = 'AssemSS';
        SSValues = zeros(13,1);
end

if (~((isw== 3||isw==6||isw==12) && numbernonlinear==0) || (exist('initializeLinKF','var') && initializeLinKF))
% Don't assemble F or K or energy if it's a linear problem

nh1 = nha;
nh2 = nhb;
nh3 = nhc;

for elem = 1:numel
    
  for ma = 1:nummat
      
   if(ieFEAP(nie-2,ma) == RegionOnElement(elem))
      
    %Extract patch material properties
    iel = MatTypeTable(2,ma); %iel   = ie(nie-1,ma); same thing;
    nonlin = MatTypeTable(3,ma);
    mateprop = MateT(ma,:);
    if iscell(mateprop) % allows for crystal plasticity or other combined material listings
        mateprop = mateprop{1};
    end
    
    if (~((isw== 3||isw==6||isw==12) && nonlin==0) || (exist('initializeLinKF','var') && initializeLinKF))
    
    %Record time of assembly
    if Compt == 1
        tic
    end
    
    
%             Compute address and offset for history variables

    ht1 = 1 + ixFEAP(1,elem) + ieFEAP(nie-3,ma);
    ht2 = 1 + ixFEAP(2,elem) + ieFEAP(nie-3,ma);
    ht3 = 1 + ixFEAP(3,elem) + ieFEAP(nie-4,ma);

%             If history variables exist move into nh1,nh2

    if(ieFEAP(nie,ma) > 0) %then
        for i = 0:ieFEAP(nie,ma)-1
          hr(nha+i) = hrvec(ht1+i);
          hr(nhb+i) = hrvec(ht2+i);
        end % i
    end

%             If Element variables exist move into nh3

    if(ieFEAP(nie-5,ma) > 0) %then
        for i = 0:ieFEAP(nie-5,ma)-1
          hr(nhc+i) = hrvec(ht3+i);
        end % i
    end
    
    %Determine element size parameters
    nel = nnz(NodesOnElement(elem,1:nen));
    nelP = getnelP(nel,ndm,nelP3,nelP4,nelP6,nelP9);
    nst = nen*ndf;
    
    %Extract patch nodal coordinates
    
    ElemFlag = NodesOnElement(elem,1:nen);
    actnode = find(ElemFlag>0);
    xl = zeros(ndm,nen);
    xl(1:ndm,actnode) = Coordinates(ElemFlag(actnode),1:ndm)';

    %Compute and Assemble Patch Stiffness
%     EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
%     [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = LocToGlobDOF2(ElemFlag, NDOFT, nel, ndf, neq);
    [EGDOFTa, EGDOFTi, ELDOFTa, ELDOFTi] = plocal(NDOFT, ElemFlag, squeeze(iedof(:,:,ma)), actnode, nen, ndf, neq);
    
    %Extract patch solution values
    ul = zeros(ndf, nen);
    ul_n = zeros(ndf, nen);
%     uld = zeros(ndf, nen);
    if transient > 0
        vl = ul;
        vl_n = ul_n;
        al = ul;
        al_n = ul_n;
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld = ul - ul_n;
%         uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
%         uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
        vl(ELDOFTa) = ModelVx(EGDOFTa)';
        vl_n(ELDOFTa) = ModelVxn_1(EGDOFTa)';
        al(ELDOFTa) = ModelAx(EGDOFTa)';
        al_n(ELDOFTa) = ModelAxn_1(EGDOFTa)';
    else
        ul(ELDOFTa) = ModelDx(EGDOFTa)';
        ul(ELDOFTi) = gBC(EGDOFTi)';
        ul_n(ELDOFTa) = ModelDxn_1(EGDOFTa)';
        ul_n(ELDOFTi) = gBC_n(EGDOFTi)';
        uld = ul - ul_n;
%         uld(ELDOFTa) = s_del_ModelDx(EGDOFTa)';
%         uld(ELDOFTi) = (gBC(EGDOFTi) - gBC_n(EGDOFTi))';
    end
    
    ElemRout
    if failFEelem
        return
    end
  
    %Assemble Element contribution to Model Quantity

%     run(AssemQuant)
    evalin('caller',[AssemQuant ';']);
    if (isw == 3 || isw == 21) && initializeLinKF && nonlin == 0
        AssemStifLL
    end
    if (isw == 5) && initializeLinKF && nonlin == 0
        AssemMassLL
    end
    
    if Compt == 1
        t = toc;
        fprintf('Element %i time to assemble: %3.4f.\n',elem,t)
    end

%             Position update terms 'ht1,ht2' from 'nh1,nh2' to save

    if(hflgu && ieFEAP(nie,ma) > 0) %then
      for i = 0:ieFEAP(nie,ma)-1
        temp      = hrvec(ht1+i);
        hrvec(ht1+i) = hr(nha+i);
        hr(nha+i) = temp;
        temp      = hrvec(ht2+i);
        hrvec(ht2+i) = hr(nhb+i);
        hr(nhb+i) = temp;
      end % i
    end

%             Position update terms 'ht3' from 'nh3' to save

    if(h3flgu && ieFEAP(nie-5,ma) > 0) %then
      for i = 0:ieFEAP(nie-5,ma)-1
        hrvec(ht3+i) = hr(nhc+i);
      end % i
    end
    
    end % if nonlin
    
   end %if ma
    
  end % ma
    
end % elem
    
end

hflgu = 0;
h3flgu = 0;
