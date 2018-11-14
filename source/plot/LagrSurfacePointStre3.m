function [S_u_v,C_u_v] = LagrSurfacePointStre3(P,C,nel,r,s,t,dim,mateprop,compon,nelP,iel)

%     -------------------------------------------------------------------

%   Compute Surface Point
    
    ind = [1 2 3 1 2 3
           1 2 3 2 3 1];
    ind1 = ind(1,compon);
    ind2 = ind(2,compon);
    Eye = eye(3);

    if nel == 4 || nel == 10
        [shp,shld,shls,be] = shltt([r,s,t],nel,nel,0,0);
        shg = shgtt(P,nel,shld,shls,nel,0,0,be);
        shpP = shltt([r,s,t],nelP,nel,0,0);
    elseif nel == 8 || nel == 27
        [shp,shld,shls,be] = shlb([r,s,t],nel,nel,0,0);
        shg = shgb(P,nel,shld,shls,nel,0,0,be);
        shpP = shlb([r,s,t],nelP,nel,0,0);
    end
    
%     for dir = 1:dim
%         for l = 1:nel
%             S_u_v(dir) = S_u_v(dir) + shp(l)*P(l,dir);
%         end
%     end
    S_u_v = P*shp;
%     for l = 1:nel
%         C_u_v = C_u_v + shp2(l)*C(l);
%     end
    C_u_v = 0;
    
    if iel == 3
        ElemE = mateprop(1);
        Elemv = mateprop(2);
        mu = ElemE/(2*(1+Elemv));
        
    for l = 1:nelP
        C_u_v = C_u_v + shpP(l)*C(4,l)*Eye(ind1,ind2) + mu*(shg(l,ind1)*C(ind2,l) + shg(l,ind2)*C(ind1,l));
    end
    for l = nelP+1:nel
        C_u_v = C_u_v + mu*(shg(l,ind1)*C(ind2,l) + shg(l,ind2)*C(ind1,l));
    end
    
    elseif iel == 1
        ElemE = mateprop(1);
        Elemv = mateprop(2);
        lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
        mu = ElemE/(2*(1+Elemv));
        
    for l = 1:nel
        C_u_v = C_u_v + lamda*(shg(l,1)*C(1,l)+shg(l,2)*C(2,l)+shg(l,3)*C(3,l))*Eye(ind1,ind2) + mu*(shg(l,ind1)*C(ind2,l) + shg(l,ind2)*C(ind1,l));
    end
    
    end