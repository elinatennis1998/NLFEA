function [tau,intb] = TauS3_PCW(xlint,xl,ul,mateprop,nel,lintt,lintl,lam)
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
    
    ib = 0;
    bf = 0;
    der = 0;
lint = lintt*lintl;
        
for l = 1:lint

    % Evaluate integ weight, determinant, and Jacobian over truncated
    % sector
    [w,ss] = intpntw(l,lintt,lintl,ib);
    [shl,shld,shls,be] = shlw(ss,6,6,der,bf);
    [QXY, shgs, Jdet, be, xs] = shgw(xl,6,shld,shls,6,bf,der,be);
        
    % Find physical location of integration point
    xy = xlint*shl;
    POUxi = POU_Coord3(xy(1),xy(2),xy(3),xl,0,nel);
    % Shape functions in physical element in order to compute displacement 
    % and FiI at integration point
    [shl,shld,shls,be] = shlb(POUxi,nel,nel,0,0);
    QXY = shgb(xl,nel,shld,shls,nel,0,0,be);
        
    % Evaluate bubble function over parent sector
    [b,db] = facebubbleW(ss);
    %Map derivatives forward to truncated sector reference configuration
    db = (db'/xs);        
%    Map forward the derivatives to spatial configuration, defined by
%    displacement field, and compute FiI
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
