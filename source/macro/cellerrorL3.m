
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
%           p = zero;
%             px = zero;
%             py = zero;
%             ux_xx = zero;
%             ux_yy = zero;
%             ux_xy = zero;
%             uy_xx = zero;
%             uy_yy = zero;
%             uy_xy = zero;
% 
%             for j = 1:nele
%                p = p + she(3,j)*ue(3,j);
%                px = px + she(1,j)*ue(3,j);
%                py = py + she(2,j)*ue(3,j);
%                ux_xx = ux_xx + shes(1,j)*ue(1,j);
%                ux_yy = ux_yy + shes(2,j)*ue(1,j);
%                ux_xy = ux_xy + shes(3,j)*ue(1,j);
%                uy_xx = uy_xx + shes(1,j)*ue(2,j);
%                uy_yy = uy_yy + shes(2,j)*ue(2,j);
%                uy_xy = uy_xy + shes(3,j)*ue(2,j);
%             end
            p = ue(3,:)*shel;
            px = ue(3,:)*she(:,1);
            py = ue(3,:)*she(:,2);
            ux_xx = ue(1,:)*shes(:,1);
            ux_yy = ue(1,:)*shes(:,2);
            ux_xy = ue(1,:)*shes(:,3);
            uy_xx = ue(2,:)*shes(:,1);
            uy_yy = ue(2,:)*shes(:,2);
            uy_xy = ue(2,:)*shes(:,3);
                            
            rx = fx + px + mu*(two*ux_xx + ux_yy + uy_xy);
            ry = fy + py + mu*(uy_xx + ux_xy + two*uy_yy);

          for k=1:nele
            for j=1:ndf

	          uf(j)  = uf(j)  + shl(k)*uc(j,k);
	          dux(j) = dux(j) + shg(k,1)*uc(j,k);
	          duy(j) = duy(j) + shg(k,2)*uc(j,k);

            end
          end
          
          uf(1) = uf(1) + (t11*rx+t12*ry)*bubble(3);
          uf(2) = uf(2) + (t21*rx+t22*ry)*bubble(3);
          dux(1) = dux(1) + (t11*rx+t12*ry)*bubble(1);
          dux(2) = dux(2) + (t21*rx+t22*ry)*bubble(1);
          duy(1) = duy(1) + (t11*rx+t12*ry)*bubble(2);
          duy(2) = duy(2) + (t21*rx+t22*ry)*bubble(2);
          

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