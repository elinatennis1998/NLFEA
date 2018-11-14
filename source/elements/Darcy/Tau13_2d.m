function [t11,t12,t21,t22] = Tau13_2d(xl,Vis,kappa,nel,nen,lint)



 
% 	t11=0;
% 	t12=0;
% 	t21=0;
% 	t22=0;
	gamma=0;
% 	beta=0;
	alpha11=0;
	alpha12=0;
	alpha21=0;
	alpha22=0;
    
    ib = 0;
    bf = 1;
    der = 0;

    for l = 1:lint
        
%         sbasis = shps(j,:);
%         s = slist(j);
%         
%         for i = 1:rnum
% 
%             w = rWgts(i)*sWgts(j);
%             rbasis = shpr(i,:);
%             r = rlist(i);
            if nel == 3 || nel == 6
                [w,litr,lits] = intpntt(l,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, xsj, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [w,litr,lits] = intpntq(l,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, xsj, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
            b = be(3);
            bx = be(1);
            by = be(2);
   
            alpha11=alpha11+Vis*(2*bx*bx+by*by)*w*xsj+Vis/kappa*b^2*w*xsj;
            alpha12=alpha12+Vis*(bx*by)*w*xsj;
            alpha21=alpha12;
            alpha22=alpha22+Vis*(bx*bx+2*by*by)*w*xsj+Vis/kappa*b^2*w*xsj;

            gamma=gamma+(b*w*xsj);

%         end %i
    end %j  ! End of Integration Loop

	beta=gamma/(alpha11*alpha22-alpha12*alpha21);

	t11= alpha22*beta;		% For stiffness terms
	t12=-alpha12*beta;
	t21=-alpha21*beta;
	t22= alpha11*beta;

