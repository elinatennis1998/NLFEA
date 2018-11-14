% 11/15/2013
% Tim Truster
%
% Array initialization for FormFE

hflgu = 0;
h3flgu = 0;

switch isw%Task Switch
    
    case 1 %Get Material Properties
        
    case 3 %Get Stiffness, Force
        hflgu = 1;
        h3flgu = 1;
        switch nreg
            case 0
                if initializeLinKF && numberlinear > 0 %initialize linear stiffness
                    KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
                    KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
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
            case 1
%                 if initializeLinKF && numberlinear > 0 %initialize linear stiffness
%                     KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
%                     KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
%                     Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
%                     Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
%                 elseif numberlinear > 0 %load linear stiffness
%                     Kdd11 = KddLL;
%                     Kdf1 = KdfLL;
%                 else %no linear stiffness, so recompute each time
                    Kdd11 = sparse(Ksdd1(:,1),Ksdd1(:,2),0,neq1,neq1);
                    Kdf1 = sparse(Ksdf1(:,1),Ksdf1(:,2),0,neq1,nieq1);
%                 end
                Kfd1 = sparse(Ksdf1(:,2),Ksdf1(:,1),0,nieq1,neq1);
                Kff1 = sparse(Ksff1(:,1),Ksff1(:,2),0,nieq1,nieq1);
                AssemQuant = 'AssemStifForc1';
            case 2
%                 if initializeLinKF && numberlinear > 0 %initialize linear stiffness
%                     KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
%                     KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
%                     Kdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
%                     Kdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
%                 elseif numberlinear > 0 %load linear stiffness
%                     Kdd11 = KddLL;
%                     Kdf1 = KdfLL;
%                 else %no linear stiffness, so recompute each time
                    Kdd22 = sparse(Ksdd2(:,1),Ksdd2(:,2),0,neq2,neq2);
                    Kdf2 = sparse(Ksdf2(:,1),Ksdf2(:,2),0,neq2,nieq2);
%                 end
                Kfd2 = sparse(Ksdf2(:,2),Ksdf2(:,1),0,nieq2,neq2);
                Kff2 = sparse(Ksff2(:,1),Ksff2(:,2),0,nieq2,nieq2);
                AssemQuant = 'AssemStifForc2';
        end
        switch transient
            case {-1,0}
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Kdd11*ModelDx - Kdf1*gBC;
                end
                Fd3 = zeros(nieq,1);
            case {1,2}
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*ModelAx - coeff1*Kdd11*ModelDx - coeff1*Kdf1*gBC;
                end
%                 F1n = F1n - alpha*Kdd11*ModelDx - alpha*Kdf1*gBC;
            case 3
                Fd1 = coeff1*Fext1;
            case 4
                Fd1 = coeff1*Fext1;
            case 5
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
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
            case 9
                switch nreg
                    case 1
                        Fd1 = Fext1;
                        if numberlinear > 0 && ~initializeLinKF
                            Fd1 = Fd1 - Kdd11*ModelDx1 - Kdf1*gBC1;
                        end
                        FdR1 = zeros(nieq1,1);
                    case 2
                        Fd2 = Fext2;
                        if numberlinear > 0 && ~initializeLinKF
                            Fd2 = Fd2 - Kdd22*ModelDx2 - Kdf2*gBC2;
                        end
                        FdR2 = zeros(nieq2,1);
                end
        end
%             if numD > 0
%                 iswtemp = isw;
%                 isw = -4;
%                 FormDN
%                 isw = iswtemp;
%             end
    case 5 %Get Mass
        Mdd11 = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
        Mdf1 = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
        Mfd = sparse(Ksdf(:,2),Ksdf(:,1),0,nieq,neq);
        Mff = sparse(Ksff(:,1),Ksff(:,2),0,nieq,nieq);
        AssemQuant = 'AssemMass';
    case 6 %Get Force
        hflgu = 1;
        h3flgu = 1;
        switch transient
            case {-1,0}
                Fd1 = Fext1;
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Kdd11*ModelDx - Kdf1*gBC;
                end
                Fd3 = zeros(nieq,1);
        AssemQuant = 'AssemForc';
            case {1,2}
%                 Fd1 = coeff1*Fext1 - F1n_1;
%                 F1n = alpha*Fext1;
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*ModelAx - coeff1*Kdd11*ModelDx - coeff1*Kdf1*gBC;
                end
%                 F1n = F1n - alpha*Kdd11*ModelDx - alpha*Kdf1*gBC;
        AssemQuant = 'AssemForc';
            case 3
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
        AssemQuant = 'AssemForc';
            case 4
                Fd1 = coeff1*Fext1;
                Fd3 = zeros(nieq,1);
        AssemQuant = 'AssemForc';
            case 5
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
                end
        AssemQuant = 'AssemForc';
            case 6
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
                end
        AssemQuant = 'AssemForc';
            case 8
                Fd1 = Fext1;
                Fd3 = zeros(nieq,1);
                if numberlinear > 0 && ~initializeLinKF
                    Fd1 = Fd1 - Mdd11*((1-Nalpham)*ModelAx+Nalpham*ModelAxn_1) - Kdd11*((1-Nalphaf)*ModelDx+Nalphaf*ModelDxn_1) - Kdf1*gBC;
                end
        AssemQuant = 'AssemForc';
            case 9
                switch nreg
                    case 1
                        Fd1 = Fext1;
                        if numberlinear > 0 && ~initializeLinKF
                            Fd1 = Fd1 - Kdd11*ModelDx1 - Kdf1*gBC1;
                        end
                        FdR1 = zeros(nieq1,1);
                        AssemQuant = 'AssemForc1';
                    case 2
                        Fd2 = Fext2;
                        if numberlinear > 0 && ~initializeLinKF
                            Fd2 = Fd2 - Kdd22*ModelDx2 - Kdf2*gBC2;
                        end
                        FdR2 = zeros(nieq2,1);
                        AssemQuant = 'AssemForc2';
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
    case 21 %Get Stiffness
        hflgu = 1;
        h3flgu = 1;
        if initializeLinKF && numberlinear > 0 %initialize linear stiffness
            KddLL = sparse(Ksdd(:,1),Ksdd(:,2),0,neq,neq);
            KdfLL = sparse(Ksdf(:,1),Ksdf(:,2),0,neq,nieq);
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