function [ omar ] = mm10_int_stress_omar( omar, mateprop, ngp )
omar.sigma_int_edge = zeros(9,ngp);
omar.sigma_int_screw = zeros(9,ngp);
omar.back_stress = zeros(12,ngp);
sigma_int_screw = 0;
sigma_int_edge = 0;
nu = mateprop.cp_prop.nu;
G = mateprop.cp_prop.mu;
b = mateprop.cp_prop.b;
R = mateprop.cp_prop.R_omar;


for ii = 1:ngp
    %% Edge
    for i=1:12
        n = omar.basis.ni(i,:,ii);
        s = omar.basis.si(i,:,ii);
        p = omar.basis.pi(i,:,ii);
        rho = omar.rho.grad18_24(i,(1:3)+3*(ii-1));
        r1 = zeros(3,3);
        r2 = zeros(3,3);
        r3 = zeros(3,3);
        r4 = zeros(3,3);
        r5 = zeros(3,3);
        for j=1:3
            for k=1:3
                for l=1:3
                    r1(k,l) =r1(k,l) + 3*rho(j)*n(j)*s(k)*s(l);
                    r2(k,l) =r2(k,l) + rho(j)*n(j)*n(k)*n(l);
                    r3(k,l) =r3(k,l) + 4*nu*rho(j)*n(j)*p(k)*p(l);
                    r4(k,l) =r4(k,l) - rho(j)*s(j)*s(k)*n(l);
                    r5(k,l) =r5(k,l) - rho(j)*s(j)*n(k)*s(l);
                end
            end
        end
        sigma_int_edge = sigma_int_edge + r1 + r2 + r3 + r4 + r5;
    end
    
    %% Screw
    for i=13:18
        n = omar.basis.ni(i,:,ii);
        s = omar.basis.si(i,:,ii);
        p = omar.basis.pi(i,:,ii);
        rho = omar.rho.grad18_24(i,(1:3)+3*(ii-1));
        r1 = zeros(3,3); r2 = zeros(3,3); r3 = zeros(3,3); r4 = zeros(3,3);
        for j=1:3
            for k=1:3
                for l=1:3
                    r1(k,l) =r1(k,l) - rho(j)*n(j)*s(k)*p(l);
                    r2(k,l) =r2(k,l) - rho(j)*n(j)*p(k)*s(l);
                    r3(k,l) =r3(k,l) + rho(j)*p(j)*s(k)*n(l);
                    r4(k,l) =r4(k,l) + rho(j)*p(j)*n(k)*s(l);
                end
            end
        end
        sigma_int_screw = sigma_int_screw + r1 + r2 + r3 + r4;
    end
    
    %%
    omar.sigma_int_edge(:,ii) = ( G*b*(R^2) )/( 8*(1-nu) ) .* reshape(sigma_int_edge,9,1);
    omar.sigma_int_screw(:,ii) = ( G*b*(R^2) )/4 .*  reshape(sigma_int_screw,9,1);
end
    omar.sum_sigma_int = omar.sigma_int_edge + omar.sigma_int_screw;
for ii = 1:ngp
    sigma_temp = reshape(omar.sum_sigma_int(:,ii),3,3);
    for i=1:12
        n = omar.basis.ni(i,:,ii)';
        s = omar.basis.si(i,:,ii);
        
        omar.back_stress(i,ii) = -s * sigma_temp * n;
    end 
end
end

