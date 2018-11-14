%% DG implementation of large deformaiton 
%% this subroutine is used for weakly enforced boundary
% pure displacement 2D case
% Tim Truster
% 10/2011
% modified by Pinlei Chen
% 04/23/2013
% for 2D large deformation with interface in it
% have tau and delta in it
% have d_ijklmn in it 
% verified for the noncomforming mesh of patch test
% for body force problem, change the lint
% change dmatL for no p
% UIUC
if isw ~= 1
CGtoDGarrays_BD
nelLP = nel;

inter = elem - (numel - numSI-numBD);
nodeAL = SurfacesD(inter,1);
nodeBL = SurfacesD(inter,2);

end


nitvms = 1;
if nitvms == 1 %VMS parameter for the stability tensor rp
pencoeff = 10;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
end


switch isw %Task Switch
%%
    case 1
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
    case 3 %interface stiffness
       % generate the initial parameter 
        ElemKLL = zeros(nst,nst);
        ElemFL = zeros(nst,1);
        NmatL = zeros(2,nst);
        BmatL = zeros(4,nst);
        BmatL1 = zeros(3,nst);        
        BmatL2 = zeros(4,nst);
        BmatL3 = zeros(4,nst);       

        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;           
        dtol = 1e-11;
        itchat = 100;12;8; %flag for iteration number to suppress chatter
             
        nelL = nel;
        xlL = xl;
        % Determin bounds of integration segment
        InterConn2D2_wBC % InterConn2DT % 
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        if nel == 3 || nel == 6  
        lint = 10;3;2;3;
        else
        lint = 10;4;10;2;3; %10 for body force problem; 4 for other problem 
        end
        ideriv = 0;
     
%         etauL = TauEE2d(xlintL,DmatL,lintt6);
%         etauR = TauEE2d(xlintR,DmatR,lintt6);
% %         tau = tauL;
%         etau = etauL + etauR;
% %         ep = tauR/tauL;
%         ep = 10*max(muL,muR)/h;
%         eb = 0;
%         Kinv = 2*eye(2);

%Form tau matrix: compute on first iteration and store; this gives
        %quadratic convergence when tau isn't recalculated (no
        %linearization of tau required).
   iterset = 4;
   if iter  <=iterset % == 0 %
        [tauL,intb] = TauS2(xl,ul,mateprop,nel,nen,lam,roL,eL1,drdr); %[Y^(-1)]

        ebL = 0;
        intedge = 0;
        
        for ie = 1:lint            
% For separate bubble types on T and Q
           if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(ie,lint,1); 
                 ebeL = edgebubble(litr,lits,nel);  %edgebubble is for T3 element
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlt(rL,lits,nelL,nel2L,0,0);
           else
                [Wgt,litr,lits] = intpntq(ie,lint,1);
%                 rL = drdr*(litr-roL)+eL1;
%                 [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nel2L,0,0);
                 ebeL = edgebubbleQ(litr,lits,nel);
           end
                                  
                    
%            b = edgebubble(litr,0);
           if nel == 3 || nel == 6
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlt(rL,lits,nel,nel,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xl,nel,shldL,shlsL,nen,0,0,be); 
           else    
               rL = drdr*(litr-roL)+eL1;
               [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
               [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xl,nel,shldL,shlsL,nen,0,0,be); 
           end

             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);            
            C1L = drdr*Wgt*Tm3L;
            
            ebL = ebL + C1L*ebeL;
            intedge = intedge + C1L;

        end   
        
        edgeK = tauL*ebL^2 ;
        ep = pencoeff*intedge*inv(edgeK); 
        ep_list(iter+1,:,inter) = [ep(1,1),ep(1,2),ep(2,1),ep(2,2)];
%         ep_List(elem,:) = [ep(1),ep(2),ep(3),ep(4),ep(5),ep(6),ep(7),ep(8),ep(9)];
%       ep_List_R(elem,:) = [tauR(1),tauR(2),tauR(3),tauR(4),tauR(5),tauR(6),tauR(7),tauR(8),tauR(9)];         
   else

        ep = [ep_list(iterset+1,1,inter) ep_list(iterset+1,2,inter)
              ep_list(iterset+1,3,inter) ep_list(iterset+1,4,inter)];
   end
%        gamL = [0.5 0
%                 0 0.5];
%        gamR = [0.5 0
%                 0 0.5];     
%        ep = 100*eye(2);
   
        ll=0;           
     for l = 1:lint

            ll = ll + 1;
            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(ll,lint,1);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(ll,lint,1);
            end
                    
%            b = edgebubble(litr,0);
            
            rL = drdr*(litr-roL)+eL1;
           if nel == 3 || nel == 6
           [shlL,shldL,shlsL,be] = shlt(rL,lits,nel,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgt(xl,nel,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgt(xl+ul,nel,shldL,shlsL,nen,0,0,be);
           else
           [shlL,shldL,shlsL,be] = shlq(rL,lits,nelL,nelL,0,0);
           [PxyL,shgsL,JdetL,bubbleL,XsL] = shgq(xl,nel,shldL,shlsL,nen,0,0,be); 
           [PxyL,shgsL,JdetL,bubbleL,xsL] = shgq(xl+ul,nel,shldL,shlsL,nen,0,0,be);      
           end
            QxyL = PxyL;  

 % initial stress sigma and material bulk term c_ijkl
           for mm = 1:nel    
           NmatL(:,2*mm-1:2*mm) = [shlL(mm,1)     0          
                                    0        shlL(mm,1)]; 
                          
           BmatL(:,2*mm-1:2*mm) = [QxyL(mm,1) 0             
                                    0         QxyL(mm,2)         
                                    QxyL(mm,2) QxyL(mm,1)        
                                     QxyL(mm,2) -QxyL(mm,1)];
           BmatL1(:,2*mm-1:2*mm) = [QxyL(mm,1) 0                  
                                     0         QxyL(mm,2)          
                                     QxyL(mm,2) QxyL(mm,1)];
           BmatL2(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
                                 QxyL(mm,2)  0                                      
                                   0     QxyL(mm,1)                
                                   0     QxyL(mm,2) ];   
           BmatL3(:,2*mm-1:2*mm)=[QxyL(mm,1)  0                     
                                   0     QxyL(mm,1)                                      
                                   QxyL(mm,2) 0                     
                                   0     QxyL(mm,2)];     
           end 
            
    %plug in the -ul to get the inverse of deformation gradient fi as 'F' in this function             
            [fiL,JxXL,FL] = kine2d(QxyL,-ul,nel,0); %this is equivalent to ikine2d
            JxXL = 1/JxXL; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
            JdetL = JdetL/JxXL;            
 
 %% weakly enforced boundary
            if exist('igbc','var')
                if igbc== 1
                    X = xl(1,1:nel)*shlL;
                    Y = xl(2,1:nel)*shlL;
                   if gbc == 1  %parabola problem
                    theta = 2.5;
                    gdisp(1) = zeros*X;               
                    gdisp(2) = theta*(X-0.25)*(X-0.75);
                   elseif gbc ==2  %
                    gdisp(1) = zeros*X;
                    gdisp(2) = zeros*Y;
                   else
                    gdisp = zeros(2,1);
                   end
                else
                    gdisp = zeros(2,1);
                end
            else
                    gdisp = zeros(2,1); 
            end 
          
            [sigma2L, cmatL] = SigmaCmat2i(FL,JxXL,mateprop,lam);
            
%             SmatL = [sigma2L(1) 0  sigma2L(3)/2 sigma2L(3)/2
%                     0 sigma2L(2)  sigma2L(3)/2 -sigma2L(3)/2
%                     sigma2L(3)/2  sigma2L(3)/2 (sigma2L(2)+sigma2L(1))/4 (sigma2L(2)-sigma2L(1))/4
%                     sigma2L(3)/2 -sigma2L(3)/2 (sigma2L(2)-sigma2L(1))/4 (sigma2L(2)+sigma2L(1))/4];            
            
            SmatL1=[sigma2L(1), sigma2L(3)
                    sigma2L(3), sigma2L(2)];           
              

             %Evaluate tangent and normal vectors
            t1L = [xsL(:,1); 0];
            [tm1L, tu1L] = VecNormalize(t1L);
            t2L = [0; 0; 1];
            tm2L = 1;
            tu2L = t2L';
            t3L = VecCrossProd(t1L,t2L);
            [tm3L, tu3L] = VecNormalize(t3L);

             c1L = drdr*Wgt*tm3L;

             
             %Evaluate tangent and normal vectors
            T1L = [XsL(:,1); 0];
            [Tm1L, Tu1L] = VecNormalize(T1L);
            T2L = [0; 0; 1];
            Tm2L = 1;
            Tu2L = T2L';
            T3L = VecCrossProd(T1L,T2L);
            [Tm3L, Tu3L] = VecNormalize(T3L);

             C1L = drdr*Wgt*Tm3L;

                %normal vectors
                nLx = tu3L(1);
                nLy = tu3L(2);

                %tagent vectors
                tLx = tu1L(1);
                tLy = tu1L(2);
           
                nvectL1= [nLx 0   nLy   
                          0  nLy  nLx ];

                nvectL2 = [eye(3,3)*nLx zeros(3,3) eye(3,3)*nLy
                           zeros(3,3)    eye(3,3)*nLy eye(3,3)*nLx];
                     
                nvecL = [nLx; nLy];

%additional stiffness term
                              
           SmatnL=[SmatL1*nvecL zeros(ndm,1)
                    zeros(ndm,1) SmatL1*nvecL];

           cmatnL=(nvectL1*cmatL(1:3,1:3));

           cmatnBL=BmatL2'*[cmatnL zeros(2,3)
                            zeros(2,3)  cmatnL ];               


            term17L=SmatnL'*BmatL2;


            term18L=nvectL1*cmatL(1:3,1:3)*BmatL1;


            term28L=NmatL'*c1L*SmatL1*nvecL;   %average stress term


             rhspulL = reshape(ul,ndf*nel,1);
             jumpu = NmatL*rhspulL - gdisp';  %jumpu=u - gdisp         

%%          assume dealta1=delata2=0.5

             
             term5L=cmatnBL*[eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]'*BmatL1;

            
             sig8L2 = [jumpu(1) jumpu(2)]*cmatnL;

  

          sig8L3 = [sig8L2(1) 0 sig8L2(3) 0;
                    sig8L2(3) 0 sig8L2(2) 0;
                    0 sig8L2(1) 0 sig8L2(3);
                    0 sig8L2(3) 0 sig8L2(2)];
               
           term8L = BmatL2'*sig8L3*BmatL3;
                                                       
              term30L=NmatL'*ep*jumpu;


             [dmatL]=dmat2_no_p(JxXL,mateprop,lam);
     
             dmatL2 = [eye(3,3)*jumpu(1) eye(3,3)*jumpu(2)]*nvectL2*dmatL/JxXL;  %missing J here
             
             term7L = BmatL1'*dmatL2*BmatL1;

             
              ElemFL = ElemFL-(-c1L*(term17L'+term18L')*jumpu-term28L+C1L*term30L);

             ElemKLL = ElemKLL - c1L*(term17L'+term18L')*NmatL - c1L*NmatL'*(term17L+term18L);
 
             ElemKLL = ElemKLL + C1L*(NmatL'*ep*NmatL);
            
            %relate to jumpu terms additional terms
             ElemKLL = ElemKLL - c1L*(term5L+term5L'+term7L+term8L);  %                     

    end %lint        

%    SurfOrientStoE
%    ul(1:ndm,1:nel) = ul(1:ndm,ilist);
%    ElemF(1:ndf*nel) = ElemFL(ilist2);
%    ElemK(1:ndf*nel,1:ndf*nel) = ElemKLL(ilist2,ilist2);
     ElemF = ElemFL;
     ElemK = ElemKLL;
end