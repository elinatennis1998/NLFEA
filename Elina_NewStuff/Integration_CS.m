%Elina Geut
%Integration for CS 

if nelL == 3 || nelL == 6
    [Wgt,r,s] = intpntt(ie,lint,1);
elseif nelL == 4 || nelL == 9
    [Wgt,r,s] = intpntq(ie,lint,1);
end

r = drdr*(r-roL)+eL1;

if nelL == 3 || nelL == 6
    [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
    [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
elseif nelL == 4 || nelL == 9
    [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
    [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
end

rR = m*(r-eL2) + eR1;

if nelR == 3 || nelR == 6
    s = 0;
else %if nelR == 4
    s = -1;
end

if nelR == 3 || nelR == 6
    [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
    [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
elseif nelR == 4 || nelR == 9
    [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
    [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
end