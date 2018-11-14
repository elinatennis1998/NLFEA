function [SurfacesL,numSL] = getpressSS49(n1,n2,m1,m2,a,b,c,d)

% Function to specify loads at ends of simply supported beam

numSL = m1;
SurfacesL = zeros(numSL,7);

for i = 1:m1/2
    SurfacesL(i,:) = [2*(n1+1)*(i-1)+1 2*(n1+1)*i+1 n1/2*(i-1)+1 4 0 0 0];
end

nind = m1/2;
ind = b - n1 - 1;
for i = a+n1/2:n1/2:c
    nind = nind + 1;
    ind = ind + 2*(n1 + 1);
    SurfacesL(nind,:) = [ind+2*(n1+1) ind i 2 0 0 0];
end