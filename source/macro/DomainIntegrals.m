function [L2,H1,ufieff,ufl2,ufix,ufiy,ufiz] = DomainIntegrals(xs,ixs,us,ProbType,elements)
%
% Function to compute domain integral from a finite element mesh and return
% the L2 norm and H1 seminorm
    
% numnp = ProbType(1);
numel = ProbType(2);
% nummat = ProbType(3);
ndm = ProbType(4);
ndf = ProbType(5);
nen = ProbType(6);

ufl2 = zeros(ndf,1);
ufix = zeros(ndf,1);
ufiy = zeros(ndf,1);
ufiz = zeros(ndf,1);
ufieff = zeros(numel,1);
    
%....	set shape function flags
ib = 0;
der = 0;
bf = 0;


%%
%-----------------------------------------------------
% Loop over elements in domain
%-----------------------------------------------------
for i = 1:length(elements)
    
    elem = elements(i);

    nel = nnz(ixs(1:nen,elem));

    if ndm == 2
    lint = IntPoint(nel);
    elseif ndm == 3
    lint = IntPoint3(nel);
    end
    uc = zeros(ndf,nen);
    
    ixc = ixs(1:nen,elem);
    actnode = find(ixc>0);
    xc = zeros(ndm,nen);
    uc = zeros(ndf,nen);
    xc(1:ndm,actnode) = xs(1:ndm,ixc(actnode));
    uc(1:ndf,actnode) = us(1:ndf,ixc(actnode));

    %....	clear the element arrays
    ufl2el = zeros(ndf,1);
    ufixel = zeros(ndf,1);
    ufiyel = zeros(ndf,1);
    ufizel = zeros(ndf,1);

    %    -----------------------------------------------------
    %     Loop over integration points
    %    -----------------------------------------------------
    for l=1:lint

    %       Evaluate shape functions

    %.... Compute Local & Global Element Shape Functions
    if ndm == 2
        if nel == 3 || nel == 6
            [w,litr,lits] =  intpntt(l,lint,ib);
            [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
            [shg, shgs, det, be, sx] = shgt(xc,nel,shld,shls,nen,bf,der,be);
        else
            [w,litr,lits] =  intpntq(l,lint,ib);
            [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
            [shg, shgs, det, be, sx] = shgq(xc,nel,shld,shls,nen,bf,der,be);
        end
    elseif ndm == 3
        if nel == 4 || nel == 10
            [w,ss] =  int3d_t(l,lint,ib);
            [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
            [shg, shgs, det, be, sx] = shgtt(xc,nel,shld,shls,nen,bf,der,be);
        else
            [w,ss] =  intpntb(l,lint,ib);
            [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
            [shg, shgs, det, be, sx] = shgb(xc,nel,shld,shls,nen,bf,der,be);
        end
    end

    c1 = det*w;	
    %
    %....	clear the fine scale solutions
    %
    uf = uc*shl;
    dux = uc*shg(:,1);
    duy = uc*shg(:,2);
    if ndm == 3
    duz = uc*shg(:,3);
    end

    %	---------------------> Error Evaluation <---------------------

    %....	loop over nodal vector

    for j=1:ndf

    ufn   = c1 * ( uf(j)^2 );
    ufpnx = c1 * ( dux(j)^2 );
    ufpny = c1 * ( duy(j)^2 );
    if ndm == 3
    ufpnz = c1 * ( duz(j)^2 );
    end

    ufl2el(j) = ufl2el(j) + ufn;
    ufixel(j) = ufixel(j) + ufpnx;
    ufiyel(j) = ufiyel(j) + ufpny;
    if ndm == 3
    ufizel(j) = ufizel(j) + ufpnz;
    end

    end

    %    -----------------------------------------------------
    %     End loop over integration points
    %    -----------------------------------------------------
    end

    %....	add the element contribution to the global error evaluated

    for j = 1:ndf

    ufl2(j) = ufl2(j) + ufl2el(j);
    ufix(j) = ufix(j) + ufixel(j);
    ufiy(j) = ufiy(j) + ufiyel(j);
    if ndm == 3
    ufiz(j) = ufiz(j) + ufizel(j);
    end

    end

    ufieff(elem) = ufieff(elem) + ufixel(1) + ufiyel(1) + ufixel(2) + ufiyel(2);
    if ndm == 3
    ufieff(elem) = ufieff(elem) + ufixel(3) + ufiyel(3) + ufizel(1) + ufizel(2) + ufizel(3);
    end

%-----------------------------------------------------
% End loop over elements in domain
%-----------------------------------------------------
end

if ndm == 2

    ul2 = sqrt(ufl2(1)+ufl2(2));
    uh1 = sqrt(ufix(1)+ufix(2)+ufiy(1)+ufiy(2));

    %
    %....	calculate the log
    %
    dl10  = log(10.d0);
    ul2  = log(ul2) / dl10;
    uh1  = log(uh1) / dl10;
    fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2,uh1)
    
    L2 = ul2;
    H1 = uh1;

    if ndf == 3

        ul2p = sqrt(ufl2(3));
        uh1p = sqrt(ufix(3)+ufiy(3));

        %
        %....	calculate the log
        %
        dl10  = log(10.d0);
        ul2p  = log(ul2p) / dl10;
        uh1p  = log(uh1p) / dl10;
        fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2p,uh1p)
        L2(2) = ul2p;
        H1(2) = uh1p;

    end

elseif ndm == 3

    ul2 = sqrt(ufl2(1)+ufl2(2)+ufl2(3));
    uh1 = sqrt(ufix(1)+ufix(2)+ufix(3)+ufiy(1)+ufiy(2)+ufiy(3)+ufiz(1)+ufiz(2)+ufiz(3));

    %
    %....	calculate the log
    %
    dl10  = log(10.d0);
    ul2  = log(ul2) / dl10;
    uh1  = log(uh1) / dl10;
    fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2,uh1)
    
    L2 = ul2;
    H1 = uh1;

    if ndf == 4

        ul2p = sqrt(ufl2(3));
        uh1p = sqrt(ufix(3)+ufiy(3));

        %
        %....	calculate the log
        %
        dl10  = log(10.d0);
        ul2p  = log(ul2p) / dl10;
        uh1p  = log(uh1p) / dl10;
        fprintf('Implicit %i  %1.7e  %1.7e\n',numel,ul2p,uh1p)
        L2(2) = ul2p;
        H1(2) = uh1p;

    end
    
end