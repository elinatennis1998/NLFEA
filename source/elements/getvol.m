function vol = getvol(xl,nel)
% 02/29/2012
% computes tau for edge bubble over a triangle for weighted Poisson problem

vol = 0;
ib = 0;
der = 0;
bf = 0;
if nel == 3 || nel == 6
    lint = 3;
else
    lint = 4;
end

for l = 1:lint
    
    if nel == 3 || nel == 6
    [Wgt,litr,lits] = intpntt(l,lint,ib);
    [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
    [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
    else
    [Wgt,litr,lits] = intpntq(l,lint,ib);
    [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
    [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
    end
    
    vol = vol + Jdet*Wgt;
    
end
