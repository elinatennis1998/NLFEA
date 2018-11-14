% Tim Truster
% 07/08/2013
% 1-D bar element
% Tested using patchbar.m and patchframe.m
% Verified for p=1,2,3 with equal spacing; I don't think unequal spacing or
% p>3 gives the correct results; second derivatives and bubble function not
% yet verified.

PatchE = mateprop(1);
PatchA = mateprop(2);
PatchI = mateprop(3);

switch isw %Task Switch
    
    case 1
        
        if ndf > 4
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 1;
        iste = nen;
%%
    case 3

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
        ElemL =norm(diff(xl'));
        
        Kxx= PatchE*PatchA/ElemL;        % Get EA/L 
        Kyy= 12*PatchE*PatchI/(ElemL^3); % Get 12EI/(L^3)
        Kzy= 6*PatchE*PatchI/(ElemL^2);  % Get 6EI/(L^2)
        Kzz1=4*PatchE*PatchI/ElemL;      % Get 4EI/L
        Kzz2=2*PatchE*PatchI/ElemL;      % Get 2EI/L   
        
        cx = diff(xl(1,:))/ElemL;
        cy = diff(xl(2,:))/ElemL;
        cx2= (diff(xl(1,:))/ElemL)^2; % cosx^2
        cy2= (diff(xl(2,:))/ElemL)^2; % cosy^2
        cxy= diff(xl(1,:))*diff(xl(2,:))/(ElemL^2); % cosx*cosy
        % Convert from local to global coordinates
        ElemK=  [ (Kxx*cx2+Kyy*cy2), (Kxx*cxy-Kyy*cxy), -Kzy*cy, -(Kxx*cx2+Kyy*cy2),-(Kxx*cxy-Kyy*cxy), -Kzy*cy 
                  (Kxx*cxy-Kyy*cxy), (Kxx*cy2+Kyy*cx2),  Kzy*cx, -(Kxx*cxy-Kyy*cxy),-(Kxx*cy2+Kyy*cx2),  Kzy*cx
                            -Kzy*cy,            Kzy*cx,    Kzz1,             Kzy*cy,           -Kzy*cx,    Kzz2
                 -(Kxx*cx2+Kyy*cy2),-(Kxx*cxy-Kyy*cxy),  Kzy*cy,  (Kxx*cx2+Kyy*cy2), (Kxx*cxy-Kyy*cxy),  Kzy*cy 
                 -(Kxx*cxy-Kyy*cxy),-(Kxx*cy2+Kyy*cx2), -Kzy*cx,  (Kxx*cxy-Kyy*cxy), (Kxx*cy2+Kyy*cx2), -Kzy*cx
                            -Kzy*cy,            Kzy*cx,    Kzz2,            -Kzy*cy,           -Kzy*cx,    Kzz1 ];
                     
        
%         Nmat = zeros(1,ndf*nel);
%         Bmat = zeros(1,ndf*nel);
%         
%         Dmat = PatchE;
%         ulres = reshape(ul,ndf*nen,1);
% 
%         % Load Gauss Points for quadrature
%         lint = nel - 1;
% %         if nel == 2
% %             lint = lintt3;%13;
% %         elseif nel == 3
% %             lint = lintq4;
% %         elseif nel == 4
% %             lint = lintt6;%13;
% %         else
% %             lint = lintq9;
% %         end
% 
%         der = 0;
%         bf = 0;
%         ib = 0;
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
%             Bmat(1:ndf:ndf*(nel-1)+1) = shg';
%             
%             sigma = Dmat*Bmat*ulres;
%             
%             ElemF = ElemF - c1*Bmat'*sigma;
%             ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
%             
%         end %je
ElemK
        
%%
    case 6

     
            
        end %je
%%
    case 15 % body force for linear element
        
%        

%%
    case 11

        

      

%%        
    case 25 %Stress Projection2

        
        
end %Task Switch
