function [S_u_v,C_u_v] = LagrSurfacePointEbub(P,C,nel,r,s,dim)

%     -------------------------------------------------------------------

%   Compute Surface Point
    
%     S_u_v = zeros(dim,1);
%     C_u_v = 0;

    if nel == 3 || nel == 6
        shp = shlt(r,s,3,3,0,0);
        ebe = edgebubble(r,s,nel);
    elseif nel == 4 || nel == 9
        shp = shlq(r,s,4,4,0,0);
        ebe = edgebubbleQ(r,s,nel);
    end
    
    S_u_v = P'*shp;
    C_u_v = C'*ebe;
%     for dir = 1:dim
%         for l = 1:nel
%             S_u_v(dir) = S_u_v(dir) + shp(l)*P(l,dir);
%         end
%     end
%     for l = 1:nel
%         C_u_v = C_u_v + shp(l)*C(l);
%     end