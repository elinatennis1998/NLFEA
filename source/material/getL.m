function L = getL(V,D,Dy,Ddy)
% computes d_ln(be)/d_be using formulas from deSouza book, Appendix 1

% Verified 5/12/12 for biaxial and nonuniform loading against plasfd.m

% Add comments

abc = [2 3
       3 1
       1 2];
tol = 1e-10;

d = diag(D);
x = V*D*V';
xv = [x(1,1); x(2,2); x(3,3); x(1,2); x(2,3); x(3,1)];
dy = diag(Dy);
ddy = diag(Ddy);

L = zeros(6,6);
Is = diag([1 1 1 1/2 1/2 1/2]);
One = [1; 1; 1; 0; 0; 0];

E = [V(1,1)^2      V(1,2)^2      V(1,3)^2
     V(2,1)^2      V(2,2)^2      V(2,3)^2
     V(3,1)^2      V(3,2)^2      V(3,3)^2
     V(1,1)*V(2,1) V(1,2)*V(2,2) V(1,3)*V(2,3)
     V(2,1)*V(3,1) V(2,2)*V(3,2) V(2,3)*V(3,3)
     V(3,1)*V(1,1) V(3,2)*V(1,2) V(3,3)*V(1,3)];

dx2dx = [ 2*x(1,1),     0,     0,           x(1,2),             0,           x(1,3) 
              0, 2*x(2,2),     0,           x(1,2),           x(2,3),             0 
              0,     0, 2*x(3,3),             0,           x(2,3),           x(1,3) 
            x(1,2),   x(1,2),     0, x(1,1)/2 + x(2,2)/2,         x(1,3)/2,         x(2,3)/2 
              0,   x(2,3),   x(2,3),         x(1,3)/2, x(2,2)/2 + x(3,3)/2,         x(1,2)/2 
            x(1,3),     0,   x(1,3),         x(2,3)/2,         x(1,2)/2, x(1,1)/2 + x(3,3)/2];
            
if abs(d(1) - d(2))>tol && abs(d(1) - d(3))>tol && abs(d(2) - d(3))>tol
    
    for a = 1:3
        
        b = abc(a,1);
        c = abc(a,2);
        
        xa = d(a);
        xb = d(b);
        xc = d(c);
        yfac = dy(a)/((xa-xb)*(xa-xc));
        
        L = L + yfac*(dx2dx - (xb+xc)*Is - (2*xa-xb-xc)*E(:,a)*E(:,a)' ...
                      - (xb-xc)*(E(:,b)*E(:,b)'-E(:,c)*E(:,c)')) ...
              + ddy(a)*E(:,a)*E(:,a)';
    end
    
elseif abs(d(1) - d(2))<=tol && abs(d(1) - d(3))<=tol
    
    L = ddy(1)*Is;
    
elseif abs(d(2) - d(3))<=tol
    
    a = 1;
    c = abc(a,2);
    xa = d(a);
    xc = d(c);
    
    xaxc = (xa-xc);
    xaxc2 = xaxc^2;
    yayc = (dy(a)-dy(c));
    s1 = yayc/xaxc2 - ddy(c)/xaxc;
    s2 = 2*xc*yayc/xaxc2 - (xa+xc)/xaxc*ddy(c);
    s3 = 2*yayc/(xaxc2*xaxc) - (ddy(a)+ddy(c))/xaxc2;
    s4 = xc*s3;
%     s5 = s4;
    s6 = xc*s4;
    
    L = s1*dx2dx - s2*Is - s3*(xv*xv') + s4*(xv*One' + One*xv') - s6*(One*One');
    
elseif abs(d(1) - d(3))<=tol
    
    a = 2;
    c = abc(a,2);
    xa = d(a);
    xc = d(c);
    
    xaxc = (xa-xc);
    xaxc2 = xaxc^2;
    yayc = (dy(a)-dy(c));
    s1 = yayc/xaxc2 - ddy(c)/xaxc;
    s2 = 2*xc*yayc/xaxc2 - (xa+xc)/xaxc*ddy(c);
    s3 = 2*yayc/(xaxc2*xaxc) - (ddy(a)+ddy(c))/xaxc2;
    s4 = xc*s3;
%     s5 = s4;
    s6 = xc*s4;
    
    L = s1*dx2dx - s2*Is - s3*(xv*xv') + s4*(xv*One' + One*xv') - s6*(One*One');
    
elseif abs(d(1) - d(2))<=tol
    
    a = 3;
    c = abc(a,2);
    xa = d(a);
    xc = d(c);
    
    xaxc = (xa-xc);
    xaxc2 = xaxc^2;
    yayc = (dy(a)-dy(c));
    s1 = yayc/xaxc2 - ddy(c)/xaxc;
    s2 = 2*xc*yayc/xaxc2 - (xa+xc)/xaxc*ddy(c);
    s3 = 2*yayc/(xaxc2*xaxc) - (ddy(a)+ddy(c))/xaxc2;
    s4 = xc*s3;
%     s5 = s4;
    s6 = xc*s4;
    
    L = s1*dx2dx - s2*Is - s3*(xv*xv') + s4*(xv*One' + One*xv') - s6*(One*One');
    
end