
%....	clear the element arrays
            ufl2el = zeros(ndf,1);
            ufixel = zeros(ndf,1);
            ufiyel = zeros(ndf,1);

%    -----------------------------------------------------
%     Loop over integration points
%    -----------------------------------------------------
            for l=1:lint

%       Evaluate shape functions

%.... Compute Local & Global Element Shape Functions
            if (neleB==3)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgt(xc,nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det, be, sx] = shgq(xc,nel,shld,shls,nen,bf,der,be);
            end
	
            c1 = det*w;	

%.....---------------------------------------------------
%     Evaluate partition of unity, uno
%.....---------------------------------------------------
% Evaluate R,S of element in terms of r,s of cell
            bigR = MR*litr + BR;
            bigS = MS*lits + BS;

% Compute element shape functions Nbar(R,S)

            if(nel==3||nel==6)
              [shel,shed,shels,bubble] = shlt(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgt(xe,nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [she,shes,ted,bubble] = shgq(xe,nel,shed,shels,nen,bf,der,bubble);
            end
%
%....	clear the fine scale solutions
%
          uf = zeros(ndf,1);
          dux = zeros(ndf,1);
          duy = zeros(ndf,1);

          for k=1:nele
            for j=1:ndf

	          uf(j)  = uf(j)  + shl(k)*uc(j,k);
	          dux(j) = dux(j) + shg(k,1)*uc(j,k);
	          duy(j) = duy(j) + shg(k,2)*uc(j,k);

            end
          end
          

%	---------------------> Error Evaluation <---------------------

%....	loop over nodal vector

          for j=1:ndf

	      ufn   = c1 * ( uf(j)^2 );
	      ufpnx = c1 * ( dux(j)^2 );
	      ufpny = c1 * ( duy(j)^2 );

	      ufl2el(j) = ufl2el(j) + ufn;
	      ufixel(j) = ufixel(j) + ufpnx;
	      ufiyel(j) = ufiyel(j) + ufpny;

          end

%    -----------------------------------------------------
%     End loop over integration points
%    -----------------------------------------------------
            end

%....	add the element contribution to the global error evaluated

       for j = 1:ndf

	      utl2(j) = utl2(j) + ufl2el(j);
	      utix(j) = utix(j) + ufixel(j);
	      utiy(j) = utiy(j) + ufiyel(j);

       end

	   utieff(elem) = utieff(elem) + ufixel(1) + ufiyel(1) + ufixel(2) + ufiyel(2);
