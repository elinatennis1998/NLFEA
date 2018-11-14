function [tau,intb] = TauDG_B(xl,D,lint,nel)
%
% 03/05/2014
% Tim Truster
% function to compute face tau for pyaramids and bricks, linear elasticity
	intb=0;
    tau = zeros(3,3);
    
    ib = 0;
    bf = 1;
    der = 0;
        
    for l = 1:lint

        %Evaluate first derivatives of basis functions at int. point
%         if nel == 4 || nel == 10
%           [w,ss] =  int3d_t(l,lint,ib);
%           [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
%           [QXY,shgs,Jdet,be,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
% 
%         else
          [w,ss] =  intpntb(l,lint,ib);
          [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);

%         end
        
%         if nel == 4 || nel == 10
%         [b,db] = facebubbleT(ss)    
%         else
        [b,db] = facebubbleH(ss);   
%         end
        bee = (db'/xs);       
         
        xsj = Jdet;        
        Bmat = [bee(1)  0       0
                0       bee(2)  0
                0       0       bee(3)
                bee(2)  bee(1)  0
                0       bee(3)  bee(2)
                bee(3)  0       bee(1)];
                
            
            tau = tau + w*xsj*(Bmat'*D*Bmat);   % for Y 

            intb = intb + xsj*w*b;

    end %j  ! End of Integration Loop
tau = inv(tau);    
