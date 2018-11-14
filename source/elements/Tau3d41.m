function [tau,intb] = Tau3d41(xl,D,lint,nel)
% 09/24/2015
% computes tau for edge bubble for von Mises plasticity, NL_Elem53_3d
%
% uses proper treatment of truncated sectors
%
% Input is history and displacement from coarse scales
% Output is stability tensor and updated history for coarse scales

% Set Material Properties

	intb=0;
    tau = zeros(3,3);
    
    ib = 0;
    bf = 1;
    der = 0;

% Loop to compute bubble function energy integral
for ll = 1:lint
    
    % Evaluate integ weight, determinant, and Jacobian over truncated
    % sector
    if nel == 4 || nel == 10
      [w,ss] =  int3d_t(ll,lint,ib);
      [~,shld,shls,be] = shltt(ss,nel,nel,der,bf);
      [~,~,Jdet,~,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);

    else
      [w,ss] =  intpntb(ll,lint,ib);
      [~,shld,shls,be] = shlb(ss,nel,nel,der,bf);
      [~,~,Jdet,~,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);

    end
    
    % Evaluate bubble function over master sector
    if nel == 4 || nel == 10
    [b,db] = facebubbleT(ss,nel); 
    else
    [b,db] = facebubbleQ(ss,nel);   
    end
    %Map derivatives forward to truncated sector reference configuration
    bee = (db'/xs);
         
	xsj = Jdet;   
        
    Bmatb = [bee(1)  0       0
            0       bee(2)  0
            0       0       bee(3)
            bee(2)  bee(1)  0
            0       bee(3)  bee(2)
            bee(3)  0       bee(1)];
            
    tau = tau + w*xsj*(Bmatb'*D*Bmatb);   % for Y 

    intb = intb + xsj*w*b;


end

tau = inv(tau);
