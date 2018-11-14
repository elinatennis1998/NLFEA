function [tau intb] = TauS2(xlint,xl,ul,mateprop,nel,nen,lam)
% 02/29/2012
% computes tau for edge bubble for weighted Poisson problem
%
% modified 01/26/2014 for proper treatment of truncated sectors

tau = 0;
ib = 0; % area integration for 2D, VOLUME intergration for 3D
der = 0;
bf = 1;
intb = 0;
if nel == 3 % Matches with Table 3 in the paper
    lint = 4;
elseif nel == 4
    lint = 9;
elseif nel == 6
    lint = 13;
else
    lint = 16;100;25;100;25;100; %100 for body force problem 
end
for ll = 1:lint
    
    % Evaluate integ weight, determinant, and Jacobian over truncated
    % sector
    if nel == 3 || nel == 6
    [w,litr,lits] = intpntt(ll,lint,ib);
    [shl,shld,shls,be] = shlt(litr,lits,3,3,der,bf);
    [QXY, shgs, Jdet, be, xs] = shgt(xlint,3,shld,shls,3,bf,der,be);
    elseif nel == 4 || nel == 9
    [w,litr,lits] = intpntq(ll,lint,ib);
    [shl,shld,shls,be] = shlq(litr,lits,4,4,der,bf);
    [QXY, shgs, Jdet, be, xs] = shgq(xlint,4,shld,shls,4,bf,der,be); 
    end
    
    % Find physical location of integration point
    xy = xlint*shl;
    POUxi = POU_Coord(xy(1),xy(2),xl,0,nel);
    % Shape functions in physical element in order to compute displacement 
    % and FiI at integration point
    if nel == 3 || nel == 6
    [shl,shld,shls,be] = shlt(POUxi(1),POUxi(2),nel,nel,der,bf);
    QXY = shgt(xl,nel,shld,shls,nen,bf,der,be);
    elseif nel == 4 || nel == 9
    [shl,shld,shls,be] = shlq(POUxi(1),POUxi(2),nel,nel,der,bf);
    QXY = shgq(xl,nel,shld,shls,nen,bf,der,be); 
    end
    
    % Evaluate bubble function over parent sector
    if nel == 3 || nel == 6
    [b,db] = edgebubble(litr,lits,nel);
    else
    [b,db] = edgebubbleQ(litr,lits,nel);
    end
    %Map derivatives forward to truncated sector reference configuration
    db = (db'/xs);
%     bee = (db'/xs)';
%     bee = db;     

%    Map forward the derivatives to spatial configuration, defined by
%    displacement field, and compute FiI
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
                
            
       [sigma2, cmat] = SigmaCmat2(F,JxX,mateprop,lam);
            
       Smat = [sigma2(1) 0  sigma2(3)/2 sigma2(3)/2
                 0 sigma2(2)  sigma2(3)/2 -sigma2(3)/2
                 sigma2(3)/2  sigma2(3)/2 (sigma2(2)+sigma2(1))/4 (sigma2(2)-sigma2(1))/4
                 sigma2(3)/2 -sigma2(3)/2 (sigma2(2)-sigma2(1))/4 (sigma2(2)+sigma2(1))/4];    
            
      Smat = Smat + cmat;
            
      tau = tau + w*xsj*(Bmat'*Smat*Bmat);   % for Y 

      intb = intb + xsj*w*b;


end

tau = inv(tau);
