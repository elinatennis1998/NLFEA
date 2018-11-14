function [t11,t12,t13,t21,t22,t23,t31,t32,t33] = TauS3(xl,ul,mateprop,nel,nen,nelP,lint)
%
% 01/04/2012
% Master tau routine for symmetric mixed form, 3D
%

    I1 = [1; 1; 1; 0; 0; 0; 0; 0; 0];
    spvec0 = I1;
    cpmat1 = I1*I1';
    cpmat2 = diag([-2,-2,-2,-1,-1,-1,0,0,0]);
    cpmat0 = cpmat1 + cpmat2;
% 	t11=0;
% 	t12=0;
% 	t21=0;
% 	t22=0;
	gamma=0;
% 	beta=0;
    alpha = zeros(3,3);
    
    ib = 0;
    bf = 1;
    der = 0;
        
    for l = 1:lint

        %Evaluate first derivatives of basis functions at int. point
        if nel == 4 || nel == 10
          [w,ss] =  int3d_t(l,lint,ib);
          [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
          shlp = shltt(ss,nelP,nel,0,0);
        else
          [w,ss] =  intpntb(l,lint,ib);
          [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
          [QXY,shgs,Jdet,be] = shgb(xl,nel,shld,shls,nel,bf,der,be);
          shlp = shlb(ss,nelP,nel,0,0);
        end
        
        [F,JxX,fi,Qxy,bee] = kine3d2(QXY,ul,nel,1,be(1:3),1);
        b = be(4);
        xsj = Jdet;
        [sigmai, cmati] = SigmaCmatNSCST3i(F,JxX,mateprop);
        [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
        spmat = JxX*theta1*spvec0;
        cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
        
        Bmat = [bee(1)  0       0
                0       bee(2)  0
                0       0       bee(3)
                bee(2)  bee(1)  0
                0       bee(3)  bee(2)
                bee(3)  0       bee(1)
                bee(2) -bee(1)  0
                0       bee(3) -bee(2)
               -bee(3)  0       bee(1)];
        
        press = ul(4,:)*shlp;
        sigmap = press*spmat;
        cmatp = press*cpmat;
        sigma2 = sigmai + sigmap;
        cmat = cmati + cmatp;
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
            
        
        alpha = alpha + w*xsj*Bmat'*(Smat+cmat)*Bmat;

        gamma=gamma+(b*w*xsj);

%         end %i
    end %j  ! End of Integration Loop
    
      ad(1,1) = alpha(2,2)*alpha(3,3) - alpha(2,3)*alpha(3,2);
      ad(1,2) = alpha(3,2)*alpha(1,3) - alpha(3,3)*alpha(1,2);
      ad(1,3) = alpha(1,2)*alpha(2,3) - alpha(1,3)*alpha(2,2);
      ad(2,1) = alpha(2,3)*alpha(3,1) - alpha(2,1)*alpha(3,3);
      ad(2,2) = alpha(3,3)*alpha(1,1) - alpha(3,1)*alpha(1,3);
      ad(2,3) = alpha(1,3)*alpha(2,1) - alpha(1,1)*alpha(2,3);
      ad(3,1) = alpha(2,1)*alpha(3,2) - alpha(2,2)*alpha(3,1);
      ad(3,2) = alpha(3,1)*alpha(1,2) - alpha(3,2)*alpha(1,1);
      ad(3,3) = alpha(1,1)*alpha(2,2) - alpha(1,2)*alpha(2,1);
      
      adj  = alpha(1,1)*ad(1,1) + alpha(1,2)*ad(2,1) + alpha(1,3)*ad(3,1);
      
	beta=gamma/adj;

	t11 = ad(1,1)*beta;		% For stiffness terms
	t12 = ad(1,2)*beta;
	t21 = ad(2,1)*beta;
	t22 = ad(2,2)*beta;
	t13 = ad(1,3)*beta;
	t31 = ad(3,1)*beta;
	t23 = ad(2,3)*beta;
	t32 = ad(3,2)*beta;
	t33 = ad(3,3)*beta;
