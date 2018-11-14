function [S_u_v,C_u_v] = LagrSurfacePointStre(P,C,nel,r,s,dim,mateprop,compon,nelP,iel)

%     -------------------------------------------------------------------

%   Compute Surface Point
    
    ind = [1 2 1
           1 2 2];
    ind1 = ind(1,compon);
    ind2 = ind(2,compon);
    Eye = eye(2);

    if nel == 3 || nel == 6
        [shp,shld,shls,be] = shlt(r,s,nel,nel,0,0);
        shg = shgt(P',nel,shld,shls,nel,0,0,be);
        shpP = shlt(r,s,nelP,nel,0,0);
    elseif nel == 4 || nel == 9
        [shp,shld,shls,be] = shlq(r,s,nel,nel,0,0);
        shg = shgq(P',nel,shld,shls,nel,0,0,be);
        shpP = shlq(r,s,nelP,nel,0,0);
    end
    
    S_u_v = P'*shp;
        
    C_u_v = 0;
    
    if iel == 3
        ElemE = mateprop(1);
        Elemv = mateprop(2);
        mu = ElemE/(2*(1+Elemv));
        
    for l = 1:nelP
        C_u_v = C_u_v + shpP(l)*C(l,3)*Eye(ind1,ind2) + mu*(shg(l,ind1)*C(l,ind2) + shg(l,ind2)*C(l,ind1));
    end
    for l = nelP+1:nel
        C_u_v = C_u_v + mu*(shg(l,ind1)*C(l,ind2) + shg(l,ind2)*C(l,ind1));
    end
    
    elseif iel == 1
        ElemE = mateprop(1);
        Elemv = mateprop(2);
        lamda = Elemv*ElemE/((1+Elemv)*(1-2*Elemv));
        mu = ElemE/(2*(1+Elemv));
        
    for l = 1:nel
        C_u_v = C_u_v + lamda*(shg(l,1)*C(l,1)+shg(l,2)*C(l,2))*Eye(ind1,ind2) + mu*(shg(l,ind1)*C(l,ind2) + shg(l,ind2)*C(l,ind1));
    end
    
%     elseif iel == 2
%         
%         C_u_v = mu(compon,:)*shg(1:2,:)*C;
        
    end
    
    