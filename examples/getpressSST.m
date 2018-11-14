function [SurfacesL,numSL] = getpressSST(n1,n2,m1,m2,a,b,c,d)

% Function to specify loads at ends of simply supported beam

numSL = 2*m1;
SurfacesL = zeros(numSL,7);

for i = 1:m1
    SurfacesL(i,:) = [(n1+1)*(i-1)+1 (n1+1)*i+1 2*n1*(i-1)+1 3 0 0 0];
end

nind = m1;
ind = b;
for i = a+2*n1:4*n1:c-2*n1
    nind = nind + 1;
    ind = ind + n1 + 1;
    SurfacesL(nind,:) = [ind+n1+1 ind i 1 0 0 0];
    nind = nind + 1;
    ind = ind + n1 + 1;
    SurfacesL(nind,:) = [ind+n1+1 ind i+2*n1 2 0 0 0];
end