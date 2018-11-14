function [tau intb] = TauS2_M(xl,ul,mateprop,nel,nelP,nen,lam,roL,eL1,drdr)
% 02/29/2012
% computes tau for edge bubble for weighted Poisson problem

    I1 = [1; 1; 0; 0];
    spvec0 = I1;
    cpmat1 = I1*I1';
    cpmat2 = diag([-2,-2,-1,0]);
    cpmat0 = cpmat1 + cpmat2;
    
tau = 0;
ib = 0; % area integration for 2D, VOLUME intergration for 3D
der = 0;
bf = 1;
intb = 0;
if nel == 3 || nel == 6
    lint = 3;4;
else
    lint = 9;16;100;25;100;25;100; %100 for body force problem 
end
for ll = 1:lint
    if nel == 3 || nel == 6
    [w,litr,lits] = intpntt(ll,lint,ib);
    rL = drdr*(litr-roL)+eL1;
    [shl,shld,shls,be] = shlt(rL,lits,nel,nel,der,bf);
    [QXY, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
    elseif nel == 4 || nel == 9
    [w,litr,lits] = intpntq(ll,lint,ib);
    rL = drdr*(litr-roL)+eL1;
    [shl,shld,shls,be] = shlq(rL,lits,nel,nel,der,bf);
    [QXY, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be); 
    end
    if nelP == 3 || nelP == 6
         [shlp,shld,shls] = shlt(litr,lits,nelP,nel,0,0);
    else
         [shlp,shld,shls] = shlq(litr,lits,nelP,nel,0,0);
    end 
    if nel == 3 || nel == 6
    [b,db] = edgebubble(rL,lits,nel);
    else
    [b,db] = edgebubbleQ(rL,lits,nel);
    end
    db = (db'/xs);
%     bee = (db'/xs)';
%     bee = db;     
%    [fi,JxX,F] = kine2d(QXY,-ul,nel,0); %this is equivalent to ikine2d
%     [fi,JxX,F,Qxy,bee] = kine2d2(QXY,-ul,nel,1,be(1:2),1);
     [F,JxX,fi,Qxy,bee] = kine2d2(QXY,ul,nel,1,db,1);
%     b = be(3);
%      JxX = 1/JxX; %this is equivalent to ikine2d
    %         [F,JxX,fi] = ikine2d(Qxy,ul,nel); %this is equivalent to above two lines
      Jdet = Jdet/JxX;    
         
     xsj = Jdet;   
     
        Bmat = [bee(1)  0       
                0       bee(2)  
                bee(2)  bee(1) 
                bee(2) -bee(1)];
       [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,mateprop);       
        spmat = JxX*spvec0;
        cpmat = JxX*cpmat0;       
        press = ul(3,:)*shlp;  
        sigmap = press*spmat;
        cmatp = press*cpmat;
        sigma2 = sigmai+sigmap;
        cmat = cmati + cmatp;        
%       [sigma2, cmat] = SigmaCmat2_M(F,JxX,mateprop,lam);
            
       Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                 0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                 sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                 sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];    
            
      Smat = Smat + cmat;
            
      tau = tau + w*xsj*(Bmat'*Smat*Bmat);   % for Y 

      intb = intb + xsj*w*b;


end

tau = inv(tau);
