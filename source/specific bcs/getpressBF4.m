function [SurfacesL,numSL] = getpressBF4(n1,m1)

% Function to specify loads at ends of simply supported beam

numSL = n1 + 2*m1;
SurfacesL = zeros(numSL,7);

nind = 1;

for i = 1:n1
    SurfacesL(nind,:) = [(n1+1)*m1+i (n1+1)*m1+i+1 n1*(m1-1)+i 3 0 0 0];
    nind = nind + 1;
end

for i = 1:m1
    SurfacesL(nind,:) = [(n1+1)*(i+1) (n1+1)*i n1*i 2 0 0 0];
    nind = nind + 1;
end

for i = 1:m1
    SurfacesL(nind,:) = [(n1+1)*(i-1)+1 (n1+1)*(i)+1 n1*(i-1)+1 4 0 0 0];
    nind = nind + 1;
end
