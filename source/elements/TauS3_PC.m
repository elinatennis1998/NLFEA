function [tau,intb] = TauS3_PC(xl,ul,mateprop,nel,lint,lam)
%
% 01/04/2012
% Master tau routine for symmetric mixed form, 3D
% 

%     I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
%     spvec0 = I1;
%     cpmat1 = I1*I1';
%     cpmat2 = diag([-2,-2,-2,-1,-1,-1,0,0,0]);
%     cpmat0 = cpmat1 + cpmat2;
% 	t11=0;
% 	t12=0;
% 	t21=0;
% 	t22=0;
	intb=0;
% 	beta=0;
    tau = zeros(3,3);
    
%     ib = 5; % Incorrect: found 1/24/14: tau was integrated over edge
%     instead of interior
%     bf = 1;
%     der = 0;
    ib = 0;
    bf = 0;
    der = 0;
    if nel == 4
        lint = 4;11;36;
    elseif nel == 8
        lint = 8;1000;64;27;
    elseif nel == 10
        lint = 14;36;
    else
        lint = 27;1000;
    end
        
    for l = 1:lint

        %Evaluate first derivatives of basis functions at int. point
        if nel == 4 || nel == 10
          [w,ss] =  int3d_t(l,lint,ib);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be,xs] = shgtt(xl,nel,shld,shls,nel,bf,der,be);

        else
          [w,ss] =  intpntb(l,lint,ib);
          [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);

        end
        
        if nel == 4 || nel == 10
        [b,db] = facebubbleT(ss,nel); 
        else
        [b,db] = facebubbleQ(ss,nel);   
        end
        db = (db'/xs);        
%         [F,JxX,fi,Qxy,bee] = kine3d2(QXY,ul,nel,1,be(1:3),1); %change to current configration
         [F,JxX,fi,Qxy,bee] = kine3d2(QXY,ul,nel,1,db,1); %change to current configration%         bee = (db'/xs)';
%         b = be(4);
%         JxX = 1/JxX; %this is equivalent to ikine2d
%         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
         Jdet = Jdet/JxX;
         
        xsj = Jdet;        
        Bmat = [bee(1)  0       0
                0       bee(2)  0
                0       0       bee(3)
                bee(2)  bee(1)  0
                0       bee(3)  bee(2)
                bee(3)  0       bee(1)
                bee(2) -bee(1)  0
                0       bee(3) -bee(2)
               -bee(3)  0       bee(1)];
                
            
            [sigma2, cmat] = SigmaCmat3(F,JxX,mateprop,lam);
            
            Smat = ...
[    sigma2(1),        0,        0,           sigma2(4)/2,                 0,           sigma2(6)/2,           sigma2(4)/2,                 0,          -sigma2(6)/2
         0,    sigma2(2),        0,           sigma2(4)/2,           sigma2(5)/2,                 0,          -sigma2(4)/2,           sigma2(5)/2,                 0
         0,        0,    sigma2(3),                 0,           sigma2(5)/2,           sigma2(6)/2,                 0,          -sigma2(5)/2,           sigma2(6)/2
   sigma2(4)/2,  sigma2(4)/2,        0, sigma2(1)/4 + sigma2(2)/4,           sigma2(6)/4,           sigma2(5)/4, sigma2(2)/4 - sigma2(1)/4,           sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2,  sigma2(5)/2,           sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,           sigma2(4)/4,          -sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,           sigma2(4)/4
   sigma2(6)/2,        0,  sigma2(6)/2,           sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4,           sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4
   sigma2(4)/2, -sigma2(4)/2,        0, sigma2(2)/4 - sigma2(1)/4,          -sigma2(6)/4,           sigma2(5)/4, sigma2(1)/4 + sigma2(2)/4,          -sigma2(6)/4,          -sigma2(5)/4
         0,  sigma2(5)/2, -sigma2(5)/2,           sigma2(6)/4, sigma2(3)/4 - sigma2(2)/4,          -sigma2(4)/4,          -sigma2(6)/4, sigma2(2)/4 + sigma2(3)/4,          -sigma2(4)/4
  -sigma2(6)/2,        0,  sigma2(6)/2,          -sigma2(5)/4,           sigma2(4)/4, sigma2(1)/4 - sigma2(3)/4,          -sigma2(5)/4,          -sigma2(4)/4, sigma2(1)/4 + sigma2(3)/4];
            
            Smat = Smat + cmat;
            
            tau = tau + w*xsj*(Bmat'*Smat*Bmat);   % for Y 

            intb = intb + xsj*w*b;

%         end %i
    end %j  ! End of Integration Loop
tau = inv(tau);    
%       ad(1,1) = alpha(2,2)*alpha(3,3) - alpha(2,3)*alpha(3,2);
%       ad(1,2) = alpha(3,2)*alpha(1,3) - alpha(3,3)*alpha(1,2);
%       ad(1,3) = alpha(1,2)*alpha(2,3) - alpha(1,3)*alpha(2,2);
%       ad(2,1) = alpha(2,3)*alpha(3,1) - alpha(2,1)*alpha(3,3);
%       ad(2,2) = alpha(3,3)*alpha(1,1) - alpha(3,1)*alpha(1,3);
%       ad(2,3) = alpha(1,3)*alpha(2,1) - alpha(1,1)*alpha(2,3);
%       ad(3,1) = alpha(2,1)*alpha(3,2) - alpha(2,2)*alpha(3,1);
%       ad(3,2) = alpha(3,1)*alpha(1,2) - alpha(3,2)*alpha(1,1);
%       ad(3,3) = alpha(1,1)*alpha(2,2) - alpha(1,2)*alpha(2,1);
%       
%       adj  = alpha(1,1)*ad(1,1) + alpha(1,2)*ad(2,1) + alpha(1,3)*ad(3,1);
%       
% 	beta=gamma/adj;
% 
% 	t11 = ad(1,1)*beta;		% For stiffness terms
% 	t12 = ad(1,2)*beta;
% 	t21 = ad(2,1)*beta;
% 	t22 = ad(2,2)*beta;
% 	t13 = ad(1,3)*beta;
% 	t31 = ad(3,1)*beta;
% 	t23 = ad(2,3)*beta;
% 	t32 = ad(3,2)*beta;
% 	t33 = ad(3,3)*beta;
