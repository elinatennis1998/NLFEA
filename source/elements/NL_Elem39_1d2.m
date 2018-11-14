% Tim Truster
% 07/06/2014
%
% Multiscale dynamics using either RFB or polynomial bubble
% d-form of dynamics
% Adds FS dofs as extra nodes, in order to investigate the spectrum and
% mode shapes.
%
% Revised version with tau and other FS quantities given in input file

% Set Material Properties

Emod = mateprop(1);
Aelem = mateprop(2);
rho = mateprop(3);
% Et = mateprop(4);
% eps_y = mateprop(5);
% sig_y = Emod*eps_y;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

bubfun = 2;1; % 1=polynomial, 2=RFB

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
%         nh1 = nen*ndm*3;
        
%%
    case 3 %Compute Stiffness and Residual
        
        ElemK = zeros(ndf*nel);
        ElemM = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel/2;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel/2);
        Bmat = zeros(1,ndf*nel/2);
        BBmat = zeros(1,ndf*nel/2);
        
        tau = mateprop(4);
        intb = mateprop(5);
        Mp = mateprop(6);
        Kp = mateprop(7);
        vol = mateprop(8);
        bave = intb/vol;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel/2-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel/2-1)+1) = shgs';
            
            Kfscs = c1*bave*( - Dmat*BBmat);
            Kfsfs = c1*Kp/vol; % Kp/lint; % 
            Kcsfs = Kfscs';
            Kcscs = c1*(Bmat'*Dmat*Bmat);
            Mfscs = c1*bave*(rho*Nmat);
            Mfsfs = c1*Mp/vol; % Mp/lint; % 
            Mcsfs = Mfscs';
            Mcscs = c1*rho*(Nmat'*Nmat);
%             Kfscsstar = intb*(rho/(Nbeta*tstep^2)*Nmat - Dmat*BBmat);
%             Kfsfsstar = inv(tau);
%             Kcsfsstar = Kfscsstar';
%             Kcscsstar = 1/(Nbeta*tstep^2)*Mcscs + Kcscs;
            
            % coarse scales
            ElemK(1:nel/2,1:nel/2) = ElemK(1:nel/2,1:nel/2) + Kcscs;
            ElemK(1:nel/2,ll+nel/2) = Kcsfs;
            % fine scales
            ElemK(ll+nel/2,1:nel/2) = Kfscs;
            ElemK(ll+nel/2,ll+nel/2) = Kfsfs;
            
            % coarse scales
            ElemM(1:nel/2,1:nel/2) = ElemM(1:nel/2,1:nel/2) + Mcscs;
            ElemM(1:nel/2,ll+nel/2) = Mcsfs;
            % fine scales
            ElemM(ll+nel/2,1:nel/2) = Mfscs;
            ElemM(ll+nel/2,ll+nel/2) = Mfsfs;

        end %je
        
        ElemF = -(ElemK*ulres + ElemM*alres);
        ElemK = ElemK + 1/(Nbeta*tstep^2)*ElemM;
        
    case 5 %Compute Mass
        
        ElemK = zeros(ndf*nel);
        ElemM = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel/2;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel/2);
        Bmat = zeros(1,ndf*nel/2);
        BBmat = zeros(1,ndf*nel/2);
        
        tau = mateprop(4);
        intb = mateprop(5);
        Mp = mateprop(6);
        Kp = mateprop(7);
        vol = mateprop(8);
        bave = intb/vol;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel/2-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel/2-1)+1) = shgs';
            
            Kfscs = c1*bave*( - Dmat*BBmat);
            Kfsfs = c1*Kp/vol;
            Kcsfs = Kfscs';
            Kcscs = c1*(Bmat'*Dmat*Bmat);
            Mfscs = c1*bave*(rho*Nmat);
            Mfsfs = c1*Mp/vol;
            Mcsfs = Mfscs';
            Mcscs = c1*rho*(Nmat'*Nmat);
%             Kfscsstar = intb*(rho/(Nbeta*tstep^2)*Nmat - Dmat*BBmat);
%             Kfsfsstar = inv(tau);
%             Kcsfsstar = Kfscsstar';
%             Kcscsstar = 1/(Nbeta*tstep^2)*Mcscs + Kcscs;
            
            % coarse scales
            ElemK(1:nel/2,1:nel/2) = ElemK(1:nel/2,1:nel/2) + Kcscs;
            ElemK(1:nel/2,ll+nel/2) = Kcsfs;
            % fine scales
            ElemK(ll+nel/2,1:nel/2) = Kfscs;
            ElemK(ll+nel/2,ll+nel/2) = Kfsfs;
            
            % coarse scales
            ElemM(1:nel/2,1:nel/2) = ElemM(1:nel/2,1:nel/2) + Mcscs;
            ElemM(1:nel/2,ll+nel/2) = Mcsfs;
            % fine scales
            ElemM(ll+nel/2,1:nel/2) = Mfscs;
            ElemM(ll+nel/2,ll+nel/2) = Mfsfs;

        end %je
        
        ElemM;

    case 6 %Compute Residual
        
%         ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel/2;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel/2);
        Bmat = zeros(1,ndf*nel/2);
        BBmat = zeros(1,ndf*nel/2);
        
        tau = mateprop(4);
        intb = mateprop(5);
        Mp = mateprop(6);
        Kp = mateprop(7);
        vol = mateprop(8);
        bave = intb/vol;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel/2-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel/2-1)+1) = shgs';
            
            Kfscs = c1*bave*( - Dmat*BBmat);
            Kfsfs = c1*Kp/vol;
            Kcsfs = Kfscs';
            Kcscs = c1*(Bmat'*Dmat*Bmat);
            Mfscs = c1*bave*(rho*Nmat);
            Mfsfs = c1*Mp/vol;
            Mcsfs = Mfscs';
            Mcscs = c1*rho*(Nmat'*Nmat);
%             Kfscsstar = intb*(rho/(Nbeta*tstep^2)*Nmat - Dmat*BBmat);
%             Kfsfsstar = inv(tau);
%             Kcsfsstar = Kfscsstar';
%             Kcscsstar = 1/(Nbeta*tstep^2)*Mcscs + Kcscs;
            
            % coarse scales
            ElemK(1:nel/2,1:nel/2) = ElemK(1:nel/2,1:nel/2) + Kcscs;
            ElemK(1:nel/2,ll+nel/2) = Kcsfs;
            % fine scales
            ElemK(ll+nel/2,1:nel/2) = Kfscs;
            ElemK(ll+nel/2,ll+nel/2) = Kfsfs;
            
            % coarse scales
            ElemM(1:nel/2,1:nel/2) = ElemM(1:nel/2,1:nel/2) + Mcscs;
            ElemM(1:nel/2,ll+nel/2) = Mcsfs;
            % fine scales
            ElemM(ll+nel/2,1:nel/2) = Mfscs;
            ElemM(ll+nel/2,ll+nel/2) = Mfsfs;

        end %je
        
        ElemF = -(ElemK*ulres + ElemM*alres);

    case 21 %Compute Stiffness
        
        ElemK = zeros(ndf*nel);
        ElemM = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel/2;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel/2);
        Bmat = zeros(1,ndf*nel/2);
        BBmat = zeros(1,ndf*nel/2);
        
        tau = mateprop(4);
        intb = mateprop(5);
        Mp = mateprop(6);
        Kp = mateprop(7);
        vol = mateprop(8);
        bave = intb/vol;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel/2-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel/2-1)+1) = shgs';
            
            Kfscs = c1*bave*( - Dmat*BBmat);
            Kfsfs = c1*Kp/vol;
            Kcsfs = Kfscs';
            Kcscs = c1*(Bmat'*Dmat*Bmat);
            Mfscs = c1*bave*(rho*Nmat);
            Mfsfs = c1*Mp/vol;
            Mcsfs = Mfscs';
            Mcscs = c1*rho*(Nmat'*Nmat);
%             Kfscsstar = intb*(rho/(Nbeta*tstep^2)*Nmat - Dmat*BBmat);
%             Kfsfsstar = inv(tau);
%             Kcsfsstar = Kfscsstar';
%             Kcscsstar = 1/(Nbeta*tstep^2)*Mcscs + Kcscs;
            
            % coarse scales
            ElemK(1:nel/2,1:nel/2) = ElemK(1:nel/2,1:nel/2) + Kcscs;
            ElemK(1:nel/2,ll+nel/2) = Kcsfs;
            % fine scales
            ElemK(ll+nel/2,1:nel/2) = Kfscs;
            ElemK(ll+nel/2,ll+nel/2) = Kfsfs;
            
            % coarse scales
            ElemM(1:nel/2,1:nel/2) = ElemM(1:nel/2,1:nel/2) + Mcscs;
            ElemM(1:nel/2,ll+nel/2) = Mcsfs;
            % fine scales
            ElemM(ll+nel/2,1:nel/2) = Mfscs;
            ElemM(ll+nel/2,ll+nel/2) = Mfsfs;

        end %je

    case -1 % boundary tractions
        
        ElemF = zeros(nst,1);
        
    case 51 % Volume stress/strain homogenization
        
        ElemSS = zeros(13,1);

    case 52 % Surface strain homogenization
        
        ElemSS = zeros(13,1);
%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,numstr+1);
        ElemS2 = zeros(nel,numstr);

        Dmat = Emod;
        Bmat = zeros(1,ndf*nel/2);
        
        % Load Guass Integration Points

        if nel == 2
            lint = 2;
            nint = 1;
        elseif nel == 3
            lint = 3;
            nint = 2;
        else
            lint = 3;
            nint = 2;
        end
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nel,1);

        %Stress Loop
        
        sw = int1d(nint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        
        for ll = 1:nint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            % Form B matrix
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            epsil = Bmat*ulres;
            sigma = Dmat*epsil;
            
            ElemS2(ll,1) = sigma;

        end %je
        
        % interpolate stress at nodes
        if nel == 2
            plist = [-1 1];
        elseif nel == 3
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        else
            sqr3 = sqrt(3);
            plist = [-sqr3 0 sqr3];
        end
        
        for ll = 1:lint
            
            r = plist(1,ll);
            shpS = shl1d(r,1,nint-1); % polynomial order p = nint - 1
            
%             for stres = 1:numstr
%                 
%                 sigmas = ElemS2(1:nint,stres)'*shpS;
%                 ElemS(ll,stres) = sigmas;
%                 
%             end
            ElemS(ll,1:numstr) = (shpS*ElemS2(1:nint,:))';
            
        end
        
%         %Integration Loop
%         Vol = xl(2) - xl(1);

        for i = 1:nel
        ElemS(i,numstr+1) = 1; % use this for simple average Vol; % use this for weighted average 
        end
        
    case 60
        
        numhr = 4;
        ElemI = zeros(10,numhr);
        
    case 12 % energy
        
        ElemE = 0;
        
        ElemK = zeros(ndf*nel);
        ElemM = zeros(ndf*nel);
%         ElemF = zeros(ndf*nel,1);

        %Set integration number
        lint = nel/2;
        
        ib = 0;
        bf = 1;
        der = 1;
        
        fb = 0;
        Dmat = Emod;
        Nmat = zeros(1,ndf*nel/2);
        Bmat = zeros(1,ndf*nel/2);
        BBmat = zeros(1,ndf*nel/2);
        
        tau = mateprop(4);
        intb = mateprop(5);
        Mp = mateprop(6);
        Kp = mateprop(7);
        vol = mateprop(8);
        bave = intb/vol;

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        vlres = reshape(vl,ndf*nel,1);
        
        sw = int1d(lint);
        [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel/2-1);
        
        if(exist('iprob','var')&&iprob==1&&step==stepmax) && elem == 1
            nind = 0;
        end
        
        for ll = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            Wgt = sw(2,ll);
            [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel/2,shld(ll,:),shls(ll,:),nen,bf,der,be);
            
            c1 = Wgt*Jdet*Aelem;
            
            Nmat(1:ndf:ndf*(nel/2-1)+1) = shl(ll,:);
                
            Bmat(1:ndf:ndf*(nel/2-1)+1) = shg';
            
            BBmat(1:ndf:ndf*(nel/2-1)+1) = shgs';
            
            Kfscs = c1*bave*( - Dmat*BBmat);
            Kfsfs = c1*Kp/vol;
            Kcsfs = Kfscs';
            Kcscs = c1*(Bmat'*Dmat*Bmat);
            Mfscs = c1*bave*(rho*Nmat);
            Mfsfs = c1*Mp/vol;
            Mcsfs = Mfscs';
            Mcscs = c1*rho*(Nmat'*Nmat);
%             Kfscsstar = intb*(rho/(Nbeta*tstep^2)*Nmat - Dmat*BBmat);
%             Kfsfsstar = inv(tau);
%             Kcsfsstar = Kfscsstar';
%             Kcscsstar = 1/(Nbeta*tstep^2)*Mcscs + Kcscs;
            
            % coarse scales
            ElemK(1:nel/2,1:nel/2) = ElemK(1:nel/2,1:nel/2) + Kcscs;
            ElemK(1:nel/2,ll+nel/2) = Kcsfs;
            % fine scales
            ElemK(ll+nel/2,1:nel/2) = Kfscs;
            ElemK(ll+nel/2,ll+nel/2) = Kfsfs;
            
            % coarse scales
            ElemM(1:nel/2,1:nel/2) = ElemM(1:nel/2,1:nel/2) + Mcscs;
            ElemM(1:nel/2,ll+nel/2) = Mcsfs;
            % fine scales
            ElemM(ll+nel/2,1:nel/2) = Mfscs;
            ElemM(ll+nel/2,ll+nel/2) = Mfsfs;

        end %je
        
        ElemE = 1/2*(ulres'*ElemK*ulres + vlres'*ElemM*vlres);
        
end