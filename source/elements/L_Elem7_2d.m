% Tim Truster
% 07/08/2013
% 1-D beam element
% Tested using patchbeam.m and patchframe.m
% Verified to give correct beam stiffness matrix

PatchE = mateprop(1);
PatchI = mateprop(2);

switch isw %Task Switch
    
    case 1
            
        if ndf > 3

            for i = 2:ndf
                lie(i,1) = 0;
            end

        end
        
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
%         if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
%         else
%             Tmat = [lx ly   0   0
%                    -ly lx   0   0
%                     0  0   lx  ly
%                     0  0  -ly  lx];
%         end
        
        Dmat = PatchE*PatchI;
        ulres = Tmat*reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = 3;
%         if nel == 2
%             lint = lintt3;%13;
%         elseif nel == 3
%             lint = lintq4;
%         elseif nel == 4
%             lint = lintt6;%13;
%         else
%             lint = lintq9;
%         end

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
%         if ndf == 3
        bdofs = [2 5 3 6];
%         else
%         bdofs = [1 3 2 4];
%         end
        len = Lelem;
        hle = len*0.5d0;
        
        for je = 1:lint

            Wgt = sw(2,je);
            [shpw,shpt] = shp1dh(sw(:,je),len);
                    
            c1 = Wgt*hle;
                
            Bmat(bdofs) = [shpw(2,:) shpt(2,:)];
            
            sigma = Dmat*Bmat*ulres(1:ndf*nel);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            
        end %je
ElemF = Tmat'*ElemF;
ElemK = Tmat'*ElemK*Tmat;
        
%%
    case 6

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
%         if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
%         else
%             Tmat = [lx ly   0   0
%                    -ly lx   0   0
%                     0  0   lx  ly
%                     0  0  -ly  lx];
%         end
        
        Dmat = PatchE*PatchI;
        ulres = Tmat*reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = 3;
%         if nel == 2
%             lint = lintt3;%13;
%         elseif nel == 3
%             lint = lintq4;
%         elseif nel == 4
%             lint = lintt6;%13;
%         else
%             lint = lintq9;
%         end

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
%         if ndf == 3
        bdofs = [2 5 3 6];
%         else
%         bdofs = [1 3 2 4];
%         end
        
        for je = 1:lint

            Wgt = sw(2,je);
            len = Lelem;
            [shpw,shpt] = shp1dh(sw(:,je),len);
                    
            c1 = Wgt*Jdet;
                
            Bmat(bdofs) = [shpw(2,:) shpt(2,:)];
            
            sigma = Dmat*Bmat*ulres(1:ndf*nel);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            
        end %je
ElemF = Tmat'*ElemF;

%%
    case 11

        ElemE = zeros(numEn,1);

%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        
    case 26
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        % Transformation
        d_ex = (xl(1,2)-xl(1,1));
        d_ey = (xl(2,2)-xl(2,1));
        Lelem = sqrt((d_ex)^2 + (d_ey)^2);
        lx = d_ex/Lelem;
        ly = d_ey/Lelem;
%         if ndf == 3
            Tmat = [lx ly 0   0   0 0
                   -ly lx 0   0   0 0
                    0  0  1   0   0 0
                    0  0  0  lx  ly 0
                    0  0  0 -ly  lx 0
                    0  0  0   0   0 1];
%         else
%             Tmat = [lx ly   0   0
%                    -ly lx   0   0
%                     0  0   lx  ly
%                     0  0  -ly  lx];
%         end
        
        Dmat = PatchE*PatchI;
        ulres = Tmat*reshape(ul,ndf*nen,1);

        % Load Gauss Points for quadrature
        lint = 3;
%         if nel == 2
%             lint = lintt3;%13;
%         elseif nel == 3
%             lint = lintq4;
%         elseif nel == 4
%             lint = lintt6;%13;
%         else
%             lint = lintq9;
%         end

        der = 0;
        bf = 0;
        ib = 0;
        sw = int1d(lint);
%         if ndf == 3
        bdofs = [2 5 3 6];
%         else
%         bdofs = [1 3 2 4];
%         end
        
        for je = 1:lint

            Wgt = sw(2,je);
            len = Lelem;
            [shpw,shpt] = shp1dh(sw(:,je),len);
                    
            c1 = Wgt*Jdet;
                
            Bmat(bdofs) = [shpw(2,:) shpt(2,:)];
            
            sigma = Dmat*Bmat*ulres(1:ndf*nel);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            
        end %je
        
        ElemS(1:2*ndf) = ElemF;
        
end %Task Switch
