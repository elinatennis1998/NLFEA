function [t11,t12,t21,t22] = TauSCST(xl,ul,mateprop,nel,nen,nelP,lint)
%
% 11/03/2011
% Master tau routine for symmetric mixed form
%
%Modified 8/23 for stored tau: Jacobian of transformation is removed so
%that it is not retained in the initial tau matrix computed in the first
%iteration while JxX is updated in future iterations in the element
%routine.

    I1 = [1; 1; 0; 0];
    spvec0 = I1;
    cpmat1 = I1*I1';
    cpmat2 = diag([-2,-2,-1,0]);
    cpmat0 = cpmat1 + cpmat2;
% 	t11=0;
% 	t12=0;
% 	t21=0;
% 	t22=0;
	gamma=0;
% 	beta=0;
    alpha = zeros(2,2);
    
    ib = 0;
    bf = 1;
    der = 0;

    %Evaluate first derivatives of basis functions at int. point
    [w,litr,lits] =  intpntt(1,1,ib);
    [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
    [QXY, shgs, Jdet, be, sx] = shgt(xl,nel,shld,shls,nel,bf,der,be);

    [F,JxX,fi] = kine2d(QXY,ul,nel,0);
    xs = sx\fi;
    xsj = Jdet;
    [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
    [theta,theta1,theta2,theta3] = ThetaS(JxX,mateprop);
            
    spmat = JxX*theta1*spvec0;
    cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
        
    for l = 1:lint
        
        %Evaluate first derivatives of basis functions at int. point
        [w,litr,lits] =  intpntt(l,lint,ib);
        [b,bee] = bubbleCST(litr,lits,xs);
        shlp = shltCST(litr,lits);
        
        Bmat = [bee(1)  0       
                0       bee(2)         
                bee(2)  bee(1)        
                bee(2) -bee(1)];
        press = ul(3,:)*shlp;
        sigmap = press*spmat;
        cmatp = press*cpmat;
        sigma3 = sigmai + sigmap;
        cmat = cmati + cmatp;
        Smat = [sigma3(1) 0  sigma3(3)/2 sigma3(3)/2
                0 sigma3(2)  sigma3(3)/2 -sigma3(3)/2
                sigma3(3)/2  sigma3(3)/2 (sigma3(2)+sigma3(1))/4 (sigma3(2)-sigma3(1))/4
                sigma3(3)/2 -sigma3(3)/2 (sigma3(2)-sigma3(1))/4 (sigma3(2)+sigma3(1))/4];
        
        alpha = alpha + w*xsj*Bmat'*(Smat+cmat)*Bmat;

        gamma=gamma+(b*w*xsj);

%         end %i
    end %j  ! End of Integration Loop

	beta=gamma/(alpha(1,1)*alpha(2,2)-alpha(1,2)*alpha(2,1));

	t11= alpha(2,2)*beta;		% For stiffness terms
	t12=-alpha(1,2)*beta;
	t21=-alpha(2,1)*beta;
	t22= alpha(1,1)*beta;

