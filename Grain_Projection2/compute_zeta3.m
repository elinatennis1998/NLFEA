% computes the shape function for the 

function [zetaX,zeta,NmatM,dx,dy] = compute_zeta3(xlR,xlL,ulresM,xlresM)
% dx =  (xlR(1,2)- xlL(1,1));
% dy =  (xlR(2,2) - xlL(2,1));
dx =  abs(xlR(1,2)- xlL(1,1));
dy =  abs(xlR(2,2) - xlL(2,1));

 
     if dx == xlresM(1)
            dy=0
%             eta_zeta = [0   dx
%                         dx   0];
            eta_zeta = [dx   0
                        0   dx];
            ulzeta=ulresM(1:2,:); 

            NmatM = [dx   0    0    0 
                    0    dx   0    0];
     elseif dy == xlresM(4) 
            dx=0;
            eta_zeta = [dy   0
                        0   dy];

            ulzeta=ulresM(3:4,:);
             NmatM= [0   0    dy    0 
                      0    0   0     dy]; 
     else
       fprintf('There is error in the selection of dx or dy\n');  

    end  
     
 zetaX = [dx
         dy];
 
  zeta=  eta_zeta*ulzeta;
   m=1;          
%  NmatM= [dx   0    dy    0 
%          0    dx   0     dy];  
%             
     
%   zeta= NmatM * ulresM;      
     