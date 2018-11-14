function [t11,t12,t21,t22] = TauNS(xl,ul,mateprop,nel,nen,nelP,lint)
%
% 10/16/2011
% Master tau routine for non-symmetric mixed form
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
        
    for l = 1:lint
        
        %Evaluate first derivatives of basis functions at int. point
        if nel == 3 || nel == 6
          [w,litr,lits] =  intpntt(l,lint,ib);
          [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
          [QXY, shgs, Jdet, be] = shgt(xl,nel,shld,shls,nel,bf,der,be);
        else
          [w,litr,lits] =  intpntq(l,lint,ib);
          [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
          [QXY, shgs, Jdet, be] = shgq(xl,nel,shld,shls,nel,bf,der,be);
        end

        if nelP == 3 || nelP == 6
          shlp = shlt(litr,lits,nelP,nel,0,0);
        else
          shlp = shlq(litr,lits,nelP,nel,0,0);
        end
        
        [F,JxX,fi,Qxy,bee] = kine2d2(QXY,ul,nel,1,be(1:2),1);
        b = be(3);
        xsj = Jdet;
        [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);
        spmat = JxX*spvec0;
        cpmat = JxX*cpmat0;
        
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

