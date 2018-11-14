function [G,H] = mm10_ornl_GH(cs)


if 0
    
G = zeros(12,12);
% H = zeros(12,12);

G0 = 0.10;
G1 = 0.22;
g2 = 0.30;
g3 = 0.38;
g4 = 0.16;
g5 = 0.45;
h2 = g2;0.05;
h3 = g3;0.12;
h4 = g4;0.03;
h5 = g5;0.25;

% Diagonal block entries of G
Gdiag = [G0 G1 G1
       G1 G0 G1
       G1 G1 G0];

tvecs = zeros(3,12);
% Get tangent vectors
for i = 1:12
    tvecs(:,i) = VecCrossProd(cs.ni(i,:),cs.bi(i,:))';
end

% Fill in template with vector norms w/o the scaling coefficients
for zeta = 1:12
    for xi = 1:12
        G(xi,zeta) = abs(cs.ni(xi,:)*tvecs(:,zeta));
    end
end
H = sqrt(ones(12,12) - G.^2); % sin = sqrt(1 - cos^2)

% Add the coefficients
G442 = [g4 g5 g3
        g5 g4 g3
        g3 g3 g2];
G424 = [g4 g3 g5
        g3 g2 g3
        g5 g3 g4];
G244 = [g4 g3 g3
        g3 g4 g5
        g3 g5 g2];
H442 = [h4 h5 h3
        h5 h4 h3
        h3 h3 h2];
H424 = [h4 h3 h5
        h3 h2 h3
        h5 h3 h4];
H244 = [h4 h3 h3
        h3 h4 h5
        h3 h5 h2];
    
% Multiply individual entries
G(1:3,1:3) = Gdiag;
G(4:6,4:6) = Gdiag;
G(7:9,7:9) = Gdiag;
G(10:12,10:12) = Gdiag;
G(4:6,1:3) = G(4:6,1:3).*G442;
G(1:3,4:6) = G(1:3,4:6).*G442;
G(7:9,1:3) = G(7:9,1:3).*G424;
G(1:3,7:9) = G(1:3,7:9).*G424;
G(10:12,1:3) = G(10:12,1:3).*G244;
G(1:3,10:12) = G(1:3,10:12).*G244;
G(7:9,4:6) = G(7:9,4:6).*G244;
G(4:6,7:9) = G(4:6,7:9).*G244;
G(10:12,4:6) = G(10:12,4:6).*G424;
G(4:6,10:12) = G(4:6,10:12).*G424;
G(10:12,7:9) = G(10:12,7:9).*G442;
G(7:9,10:12) = G(7:9,10:12).*G442;

H(1:3,1:3) = Gdiag;
H(4:6,4:6) = Gdiag;
H(7:9,7:9) = Gdiag;
H(10:12,10:12) = Gdiag;
H(4:6,1:3) = H(4:6,1:3).*H442;
H(1:3,4:6) = H(1:3,4:6).*H442;
H(7:9,1:3) = H(7:9,1:3).*H424;
H(1:3,7:9) = H(1:3,7:9).*H424;
H(10:12,1:3) = H(10:12,1:3).*H244;
H(1:3,10:12) = H(1:3,10:12).*H244;
H(7:9,4:6) = H(7:9,4:6).*H244;
H(4:6,7:9) = H(4:6,7:9).*H244;
H(10:12,4:6) = H(10:12,4:6).*H424;
H(4:6,10:12) = H(4:6,10:12).*H424;
H(10:12,7:9) = H(10:12,7:9).*H442;
H(7:9,10:12) = H(7:9,10:12).*H442;

else
    
    load('GHmat_ornl.mat');

end
