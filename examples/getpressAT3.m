function [SurfacesL,numSL] = getpressAT3(n1,m1,p)

% Function to specify loads at ends of simply supported beam

numSL = m1;
SurfacesL = zeros(numSL,7);

for i = 1:m1
    SurfacesL(i,:) = [(n1+1)*(i+1) (n1+1)*i 2*n1*i 1 p 0 0];
end
