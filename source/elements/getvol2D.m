function vol = getvol2D(xl,lints,nel)
% 09/11/2012

if nel == 3
    lint = lints(1);%13;
elseif nel == 4
    lint = lints(2);
elseif nel == 6
    lint = lints(3);%13;
elseif nel == 9
    lint = lints(4);
end
vol = 0;
ib = 0;
der = 0;
bf = 0;

for l = 1:lint
    
    if nel == 3 || nel == 6
        [Wgt,litr,lits] = intpntt(l,lint,ib);
        [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
        [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nel,bf,der,be);
    elseif nel == 4 || nel == 9
        [Wgt,litr,lits] = intpntq(l,lint,ib);
        [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
        [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nel,bf,der,be);
    end
    
    vol = vol + Jdet*Wgt;
    
end
