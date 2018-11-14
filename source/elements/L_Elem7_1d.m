% Tim Truster
% 07/08/2013
% 1-D beam element
% Tested using patchbeam.m and patchframe.m
% Verified to give correct beam stiffness matrix

PatchE = mateprop(1);
PatchI = mateprop(2);

switch isw %Task Switch
    
    case 1
        
        if nen > 2 % zero out interior nodes n>2, meaning n1 and n2 are the endpoints
            
            for j = 4:nen+1
                lie(i,j) = 0;
            end
            
            if ndf == 3

                for j = 2:3
                for i = 1
                    lie(i,j) = 0;
                end
                end

            elseif ndf > 3

                for j = 2:3
                for i = 2:ndf
                    lie(i,j) = 0;
                end
                end

            end
            
        else
            
            if ndf == 3

                for i = 1
                    lie(i,1) = 0;
                end

            elseif ndf > 3

                for i = 2:ndf
                    lie(i,1) = 0;
                end

            end
        
        end
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        Dmat = PatchE*PatchI;
        ulres = reshape(ul,ndf*nen,1);

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
        if ndf == 3
        bdofs = [2 5 3 6];
        else
        bdofs = [1 3 2 4];
        end
        len = xl(1,2) - xl(1,1);
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
ElemK;
        
%%
    case 6

        ElemF = zeros(ndf*nel,1);
        Nmat = zeros(1,ndf*nel);
        Bmat = zeros(1,ndf*nel);
        
        Dmat = PatchE;
        ulres = reshape(ul,ndf*nen,1);

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
        if ndf == 3
        bdofs = [2 5 3 6];
        else
        bdofs = [1 3 2 4];
        end
        
        for je = 1:lint

            Wgt = sw(2,je);
            len = xl(1,2) - xl(1,1);
            [shpw,shpt] = shp1dh(sw(:,je),len);
                    
            c1 = Wgt*Jdet;
                
            Bmat(bdofs) = [shpw(2,:) shpt(2,:)];
            
            sigma = Dmat*Bmat*ulres*ulres(1:ndf*nel);
            
            ElemF = ElemF - c1*Bmat'*sigma;
            
        end %je
%%
    case 15 % body force for linear element
        
% %         ElemF = zeros(nst,1); In FormFE
%         
%         Nmat = zeros(1,nst);
%         
%         % Load Gauss Points for quadrature
%         lint = 3;
% %         if nel == 3
% %             lint = lintt3;%13;
% %         elseif nel == 4
% %             lint = lintq4;
% %         elseif nel == 6
% %             lint = lintt6;%13;
% %         elseif nel == 9
% %             lint = lintq9;
% %         end
% 
%         der = 0;
%         bf = 0;
%         
%         if exist('iprob','var') && iprob == 1
%             fb = 10;
%         else
%             fb = bodyf;
%         end
%         
%         sw = int1d(lint);
%         [shl,shld,shls,be] = shl1d(sw(1,:),lint,nel-1);
%         
%         for je = 1:lint
% 
%             Wgt = sw(2,je);
%             [shg, shgs, Jdet, be] = shg1d(xl,ndm,nel,shld(je,:),shls(je,:),nen,bf,der,be);
%                     
%             c1 = Wgt*Jdet*PatchA;
%                 
%             Nmat(1:ndf:ndf*(nel-1)+1) = shl(je,:);
%             
%           case 25 % Stress averaging
%           ElemF = ElemF + c1*(Nmat'*fb);
%                 
%         end %je

%%
    case 11

        ElemE = zeros(numEn,1);

%%        
    case 25 %Stress Projection2

        ElemS = zeros(nel,npstr+1);
        ElemS2 = zeros(nel,npstr);
        
end %Task Switch
