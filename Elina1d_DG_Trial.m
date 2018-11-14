 %Elina Geut
 %DG subroutine, first trial
 %Created 10/20/2018


        ElemF = zeros((nelL+nelR)*ndf,1);
        
        PatchEL = MateT(mateprop(1),1);
        PatchAL = MateT(mateprop(1),2);
        PatchER = MateT(mateprop(2),1);
        PatchAR = MateT(mateprop(2),2);

        DmatL = PatchEL*PatchAL;
        DmatR = PatchER*PatchAR;
        
%         rL = 1;                                                            %Param. coord. of L.Node 
%         rR = -1;                                                           %Param. coord. of R.Node  
%         
%         nL = -1;                                                           %Unit normal of L.S. ---> (-)
%         nR = 1;                                                            %Unit normal of R.S. <--- (+)
%         
%         NparamL = [(1-rL)/2 (1+rL)/2];                                     %Shape func. in parametric form
%         NparamR = [(1-rR)/2 (1+rR)/2];                                      
%         DerBL = [-1/2 1/2];
%         DerBR = [-1/2 1/2];
%         
%         %Get Jacobain that will relate parametric to master
%         
%         JL = DerBL*xl(1:nelL)';                                            %Jacobian for L.Elem
%         JR = DerBR*xl(nelL+1:nelL+nelR)';                                  %Jacobian for R.Elem
%         
%         NmatL = NparamL;
%         NmatR = NparamR;
%         
%         BmatL = (1/JL)*DerBR;
%         BmatR = (1/JR)*DerBL;
%         
%         %ngp = (p+1)/2=1
%         
%         w = 1;                                                             %weights for gauss pts
%         DetJL = det(JL);
%         DetJR = det(JR);
% 
%        %Second Term of Weak Form 
%         
%         KLL2 = (-0.5*w).*NmatL'*DmatL*BmatL; 
%         KLR2 = (-0.5*w).*NmatL'*DmatR*BmatR;
%         KRL2 = (0.5*w).*NmatR'*DmatL*BmatL;
%         KRR2 = (0.5*w).*NmatR'*DmatR*BmatR;
%         
%         K2 = [KLL2 KLR2                                                      
%               KRL2 KRR2];
%         %Third Term of Weak Form
%         
%         K3 = K2';
%         
%          %Third Term of Weak Form       
% 
%         tau = 10^(1);  
%         KLL4 = (0.5*w*tau).*(NmatL'*NmatL);
%         KLR4 = (-0.5*w*tau).*(NmatL'*NmatR);
%         KRL4 = (-0.5*w*tau).*(NmatR'*NmatL);
%         KRR4 = (0.5*w*tau).*(NmatR'*NmatR);
%         
%         
%         K4 = [KLL4 KLR4                                                      
%               KRL4 KRR4];
%         
%        ElemK = K2+K3+K4;
%        
       
         % Use one point integration rule
        
        lintL;lintR;ndm;ndf;pL;pR;nel;