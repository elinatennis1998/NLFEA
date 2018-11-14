zero = 0.d0;
%     No modifications when copied to NLFEA ver2
%....	clear the element arrays
            ufl2el = zeros(ndf,1);
            ufixel = zeros(ndf,1);
            ufiyel = zeros(ndf,1);
        I1 = [1; 1; 0; 0];
        spvec0 = I1;
        cpmat1 = I1*I1';
        cpmat2 = diag([-2,-2,-1,0]);
        cpmat0 = cpmat1 + cpmat2;
        nel = nele;

%    -----------------------------------------------------
%     Loop over integration points
%    -----------------------------------------------------
            for l=1:lint

%       Evaluate shape functions

%.... Compute Local & Global Element Shape Functions
            if(nel==3||nel==6)
              [w,litr,lits] =  intpntt(l,lint,ib);
              [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
              [shg, shgs, det] = shgt(xc(:,1:nel),nel,shld,shls,nen,bf,der,be);
            else
              [w,litr,lits] =  intpntq(l,lint,ib);
              [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
              [shg, shgs, det] = shgq(xc(:,1:nel),nel,shld,shls,nen,bf,der,be);
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
              [Qxye,shes,ted,bubble] = shgt(xe+ue(1:2,:),nel,shed,shels,nen,bf,der,bubble);
            else
              [shel,shed,shels,bubble] = shlq(bigR,bigS,nel,nel,der,bf);
              [Qxye,shes,ted,bubble] = shgq(xe+ue(1:2,:),nel,shed,shels,nen,bf,der,bubble);
            end

            [fi,JxX,F,QXYe] = kine2d(Qxye,-ue,nel,1);
            JxX = 1/JxX;
            bubble(1:2) = bubble(1:2)*F;
            
            if symmns == 1 %NS
            spvec = JxX*spvec0;
            cpmat = JxX*cpmat0;
            press = ue(3,:)*shel;
            else %S
            [theta,theta1,theta2,theta3] = ThetaS(JxX,d);
            spvec = JxX*theta1*spvec0;
            cpmat = (theta2*JxX^2 + theta1*JxX)*cpmat1 + theta1*JxX*cpmat2;
            press = ue(3,:)*shel;
            end
            sigmap = press*spvec;
            cmatp = press*cpmat;
            [sigmai, cmati] = SigmaCmatNSCST2i(F,JxX,d);
%
%....	clear the fine scale solutions
%
          uf = zeros(ndf,1);
          dux = zeros(ndf,1);
          duy = zeros(ndf,1);
            px = zero;
            py = zero;

            for j = 1:nele
               px = px + ue(3,j)*Qxye(j,1);
               py = py + ue(3,j)*Qxye(j,2);
            end

                if iprob == 6
%                     mu = 40;
                    X = xe(1,:)*shel;
                    Y = xe(2,:)*shel;
                    fbx = - (101*40)/(101*X + 100) - (101*40*(101*X + 100))/10000;
        if nel == 3
        % Form BB matrix
        BBmat = [zeros(6,9)        
                 0         0        Qxye(1,1) 0         0        Qxye(2,1) 0         0        Qxye(3,1) 
                 0         0        Qxye(1,2) 0         0        Qxye(2,2) 0         0        Qxye(3,2)];
        elseif nel == 4
        % Form BB matrix
        BBmat = [shes(1,1) 0         0        shes(2,1) 0         0        shes(3,1) 0         0        shes(4,1) 0         0
                 shes(1,2) 0         0        shes(2,2) 0         0        shes(3,2) 0         0        shes(4,2) 0         0
                 shes(1,3) 0         0        shes(2,3) 0         0        shes(3,3) 0         0        shes(4,3) 0         0
                 0         shes(1,1) 0        0         shes(2,1) 0        0         shes(3,1) 0        0         shes(4,1) 0
                 0         shes(1,2) 0        0         shes(2,2) 0        0         shes(3,2) 0        0         shes(4,2) 0
                 0         shes(1,3) 0        0         shes(2,3) 0        0         shes(3,3) 0        0         shes(4,3) 0
                 0         0         Qxye(1,1) 0         0        Qxye(2,1) 0         0        Qxye(3,1) 0         0        Qxye(4,1) 
                 0         0         Qxye(1,2) 0         0        Qxye(2,2) 0         0        Qxye(3,2) 0         0        Qxye(4,2)];
        elseif nel == 6
        % Form BB matrix
        BBmat = [shes(1,1) 0         0        shes(2,1) 0         0        shes(3,1) 0         0        shes(4,1) 0         0        shes(5,1) 0         0        shes(6,1) 0         0
                 shes(1,2) 0         0        shes(2,2) 0         0        shes(3,2) 0         0        shes(4,2) 0         0        shes(5,2) 0         0        shes(6,2) 0         0
                 shes(1,3) 0         0        shes(2,3) 0         0        shes(3,3) 0         0        shes(4,3) 0         0        shes(5,3) 0         0        shes(6,3) 0         0
                 0         shes(1,1) 0        0         shes(2,1) 0        0         shes(3,1) 0        0         shes(4,1) 0        0         shes(5,1) 0        0         shes(6,1) 0
                 0         shes(1,2) 0        0         shes(2,2) 0        0         shes(3,2) 0        0         shes(4,2) 0        0         shes(5,2) 0        0         shes(6,2) 0
                 0         shes(1,3) 0        0         shes(2,3) 0        0         shes(3,3) 0        0         shes(4,3) 0        0         shes(5,3) 0        0         shes(6,3) 0
                 0         0         Qxye(1,1) 0         0         Qxye(2,1) 0         0         Qxye(3,1) 0         0         Qxye(4,1) 0         0         Qxye(5,1) 0         0         Qxye(6,1)
                 0         0         Qxye(1,2) 0         0         Qxye(2,2) 0         0         Qxye(3,2) 0         0         Qxye(4,2) 0         0         Qxye(5,2) 0         0         Qxye(6,2)];
        elseif nel == 9
        % Form BB matrix
        BBmat = [shes(1,1) 0         0        shes(2,1) 0         0        shes(3,1) 0         0        shes(4,1) 0         0        shes(5,1) 0         0        shes(6,1) 0         0        shes(7,1) 0         0        shes(8,1) 0         0        shes(9,1) 0         0
                 shes(1,2) 0         0        shes(2,2) 0         0        shes(3,2) 0         0        shes(4,2) 0         0        shes(5,2) 0         0        shes(6,2) 0         0        shes(7,2) 0         0        shes(8,2) 0         0        shes(9,2) 0         0
                 shes(1,3) 0         0        shes(2,3) 0         0        shes(3,3) 0         0        shes(4,3) 0         0        shes(5,3) 0         0        shes(6,3) 0         0        shes(7,3) 0         0        shes(8,3) 0         0        shes(9,3) 0         0
                 0         shes(1,1) 0        0         shes(2,1) 0        0         shes(3,1) 0        0         shes(4,1) 0        0         shes(5,1) 0        0         shes(6,1) 0        0         shes(7,1) 0        0         shes(8,1) 0        0         shes(9,1) 0
                 0         shes(1,2) 0        0         shes(2,2) 0        0         shes(3,2) 0        0         shes(4,2) 0        0         shes(5,2) 0        0         shes(6,2) 0        0         shes(7,2) 0        0         shes(8,2) 0        0         shes(9,2) 0
                 0         shes(1,3) 0        0         shes(2,3) 0        0         shes(3,3) 0        0         shes(4,3) 0        0         shes(5,3) 0        0         shes(6,3) 0        0         shes(7,3) 0        0         shes(8,3) 0        0         shes(9,3) 0
                 0         0         Qxye(1,1) 0         0         Qxye(2,1) 0         0         Qxye(3,1) 0         0         Qxye(4,1) 0         0         Qxye(5,1) 0         0         Qxye(6,1) 0         0         Qxye(7,1) 0         0         Qxye(8,1) 0         0         Qxye(9,1)
                 0         0         Qxye(1,2) 0         0         Qxye(2,2) 0         0         Qxye(3,2) 0         0         Qxye(4,2) 0         0         Qxye(5,2) 0         0         Qxye(6,2) 0         0         Qxye(7,2) 0         0         Qxye(8,2) 0         0         Qxye(9,2)];
        end
                else
                    
                end
                
        I3 = [1 0 0 0; 0 1 0 0; 0 0 2 0];
        I4 = [1 0 0 0 0 1
              0 0 1 0 1 0];
          
            if symmns == 1 %NS
                
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
            ddu = BBmat*reshape(ue,3*nel,1);
            ddu = ddu(1:6);

            Fddu = F*ue(1:2,:)*shes;
            Fddu4 = Fddu*I3;
            diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
            bigF = [F(1,1) 0 0 F(1,2) 0 0
                    0 0 F(2,1) 0 0 F(2,2)
                    F(2,1) 0 F(1,1) F(2,2) 0 F(1,2)
                    0 0 F(1,1) 0 0 F(1,2)
                    0 F(2,1) 0 0 F(2,2) 0
                    0 F(1,1) F(2,1) 0 F(1,2) F(2,2)];
            term2 = Fddu4*sigma; %page 6a
            term3 = I4*diagc*bigF*ddu; %page 6a

            %Evaluate residual of equilibrium equation
            rx = JxX*px + term2(1) + term3(1) + fbx;
            ry = JxX*py + term2(2) + term3(2) + fby;
            
            else %S
                
            sigma = sigmai + sigmap;
            cmat = cmati + cmatp;
        ddu = BBmat*reshape(ue,3*nel,1);
        ddu = ddu(1:6);

        Fddu = F*ue(1:2,:)*shes;
        Fddu4 = Fddu*I3;
        diagc = [cmat(1:3,1:3) zeros(3); zeros(3) cmat(1:3,1:3)];
        bigF = [F(1,1) 0 0 F(1,2) 0 0
                0 0 F(2,1) 0 0 F(2,2)
                F(2,1) 0 F(1,1) F(2,2) 0 F(1,2)
                0 0 F(1,1) 0 0 F(1,2)
                0 F(2,1) 0 0 F(2,2) 0
                0 F(1,1) F(2,1) 0 F(1,2) F(2,2)];
        term2 = Fddu4*sigma; %page 6a
        term3 = I4*diagc*bigF*ddu; %page 6a

        %Evaluate residual of equilibrium equation
        rx = JxX*theta1*px + term2(1) + term3(1) + fbx;
        ry = JxX*theta1*py + term2(2) + term3(2) + fby;
                
            end  

          for k=1:nele
            for j=1:ndf

	          uf(j)  = uf(j)  + uc(j,k)*shl(k);
% 	          dux(j) = dux(j) + uc(j,k)*Qxy(k,1);
% 	          duy(j) = duy(j) + uc(j,k)*Qxy(k,2);
	          dux(j) = dux(j) + uc(j,k)*shg(k,1);
	          duy(j) = duy(j) + uc(j,k)*shg(k,2);

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