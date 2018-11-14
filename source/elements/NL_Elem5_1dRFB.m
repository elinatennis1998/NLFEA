
% Tim Truster
% 2/12/2014
% 1d linear dynamic element with a RFB computed from the reaction-diffusion
% equation as a 3rd shape function.
% The element gives improved results for CFL = 0.5 but still is not stable
% below CFL=0.1

% Set Material Properties

cwave = mateprop(1);
t_on = 1;0;

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
%%    
    case 1
        
        if ndf > 1
            
            for i = 2:ndf
                lie(i,1) = 0;
            end
            
        end
        
%         nh1 = nen*ndm*3;
        
%%
    case 3 %Compute Stiffness and Residual
        
        e = sqrt(Nbeta)*cwave*tstep;
        h = xl(2) - xl(1);
%         Me = ....
% [                                  h/3,                                  h/6,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%                                    h/6,                                  h/3,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%   (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*(3*e + 8*exp(2/e) + 2*exp(4/e) - 3*e*exp(4/e) + 2))/(2*(exp(2/e) + 1)^2)];
% Ke = ...
% [  cwave^2/h, -cwave^2/h,                                                                     0
%   -cwave^2/h,  cwave^2/h,                                                                     0
%            0,          0, -(2*cwave^2*(4*exp(2/e) - e*(exp(4/e) - 1)))/(e^2*h*(exp(2/e) + 1)^2)];
        Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
Ke = ...
[  cwave^2/h, -cwave^2/h,                                                                      0 
  -cwave^2/h,  cwave^2/h,                                                                      0 
           0,          0, (cwave^2*(e*(exp((2*h)/e) - 1) - 2*h*exp(h/e)))/(e^2*(exp(h/e) + 1)^2)];

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        ElemF = - Me*alres - Ke*ulres;
        ElemK = Me + (Nbeta*tstep^2)*Ke;
            
        if(step==stepmax) && elem == 1
            nind = 0;
        end
        if(step==stepmax)
        sw = int1d(2);
        [shl,shld,shls,be] = shl1d(sw(1,:),2,1);
        Dmat = cwave^2;
        
        for ll = 1:2                    

            % Evaluate 1-D basis functions at integration points
            [shg, shgs, Jdet, be] = shg1d(xl(1:2),ndm,2,shld(ll,:),shls(ll,:),2,0,0,0);
            
            xint = xl(1:2)*shl(ll,:)';
            a = 1/e;
            c1 = -1/(exp(a) + exp(-a));
            c2 = c1;
            du = 2*a/h*(c1*exp(a*sw(1,ll)) - c2*exp(-a*sw(1,ll)));
            Bmat(1:ndf:ndf*(nel-1)+1) = [shg' du];

            sigma = Dmat*(Bmat*ulres);
            
            nind = nind + 1;
            Sig1(nind,:)=[xint/L sigma/NodeLoadnp(3)];
                

        end %je
        end
        
        
    case 5 %Compute Mass
        
        e = sqrt(Nbeta)*cwave*tstep;
        h = xl(2) - xl(1);
%         Me = ....
% [                                  h/3,                                  h/6,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%                                    h/6,                                  h/3,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%   (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*(3*e + 8*exp(2/e) + 2*exp(4/e) - 3*e*exp(4/e) + 2))/(2*(exp(2/e) + 1)^2)];
        Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
        ElemM = Me;

    case 6 %Compute Residual
        
        e = sqrt(Nbeta)*cwave*tstep;
        h = xl(2) - xl(1);
%         Me = ....
% [                                  h/3,                                  h/6,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%                                    h/6,                                  h/3,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%   (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*(3*e + 8*exp(2/e) + 2*exp(4/e) - 3*e*exp(4/e) + 2))/(2*(exp(2/e) + 1)^2)];
% Ke = ...
% [  cwave^2/h, -cwave^2/h,                                                                     0
%   -cwave^2/h,  cwave^2/h,                                                                     0
%            0,          0, -(2*cwave^2*(4*exp(2/e) - e*(exp(4/e) - 1)))/(e^2*h*(exp(2/e) + 1)^2)];
        Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
Ke = ...
[  cwave^2/h, -cwave^2/h,                                                                      0 
  -cwave^2/h,  cwave^2/h,                                                                      0 
           0,          0, (cwave^2*(e*(exp((2*h)/e) - 1) - 2*h*exp(h/e)))/(e^2*(exp(h/e) + 1)^2)];

        % Compute kinematic fields at intermediate time levels
        ul_n_af = (1 - Nalphaf)*ul + Nalphaf*ul_n;
        al_n_am = (1 - Nalpham)*al + Nalpham*al_n;
        ulres = reshape(ul_n_af,ndf*nel,1);
        alres = reshape(al_n_am,ndf*nel,1);
        
        ElemF = - Me*alres - Ke*ulres;

    case 21 %Compute Stiffness
        
        e = sqrt(Nbeta)*cwave*tstep;
        h = xl(2) - xl(1);
%         Me = ....
% [                                  h/3,                                  h/6,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%                                    h/6,                                  h/3,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%   (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*(3*e + 8*exp(2/e) + 2*exp(4/e) - 3*e*exp(4/e) + 2))/(2*(exp(2/e) + 1)^2)];
% Ke = ...
% [  cwave^2/h, -cwave^2/h,                                                                     0
%   -cwave^2/h,  cwave^2/h,                                                                     0
%            0,          0, -(2*cwave^2*(4*exp(2/e) - e*(exp(4/e) - 1)))/(e^2*h*(exp(2/e) + 1)^2)];
        Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
Ke = ...
[  cwave^2/h, -cwave^2/h,                                                                      0 
  -cwave^2/h,  cwave^2/h,                                                                      0 
           0,          0, (cwave^2*(e*(exp((2*h)/e) - 1) - 2*h*exp(h/e)))/(e^2*(exp(h/e) + 1)^2)];
        
        ElemK = Me + (Nbeta*tstep^2)*Ke;
        
    case 40 % Initialize FS acceleration
        
        
    case 12 % energy
        
        e = sqrt(Nbeta)*cwave*tstep;
        h = xl(2) - xl(1);
%         Me = ....
% [                                  h/3,                                  h/6,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%                                    h/6,                                  h/3,                                        (h*((2*e)/(exp(2/e) + 1) - e + 1))/2
%   (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*((2*e)/(exp(2/e) + 1) - e + 1))/2, (h*(3*e + 8*exp(2/e) + 2*exp(4/e) - 3*e*exp(4/e) + 2))/(2*(exp(2/e) + 1)^2)];
% Ke = ...
% [  cwave^2/h, -cwave^2/h,                                                                     0
%   -cwave^2/h,  cwave^2/h,                                                                     0
%            0,          0, -(2*cwave^2*(4*exp(2/e) - e*(exp(4/e) - 1)))/(e^2*h*(exp(2/e) + 1)^2)];
        Me = ...
[                            h/3,                            h/6,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
                             h/6,                            h/3,                                                h/2 - e + (2*e)/(exp(h/e) + 1) 
  h/2 - e + (2*e)/(exp(h/e) + 1), h/2 - e + (2*e)/(exp(h/e) + 1), (3*e + h - 3*e*exp((2*h)/e) + 4*h*exp(h/e) + h*exp((2*h)/e))/(exp(h/e) + 1)^2];
Ke = ...
[  cwave^2/h, -cwave^2/h,                                                                      0 
  -cwave^2/h,  cwave^2/h,                                                                      0 
           0,          0, (cwave^2*(e*(exp((2*h)/e) - 1) - 2*h*exp(h/e)))/(e^2*(exp(h/e) + 1)^2)];

        ulres = reshape(ul,ndf*nel,1);
        vlres = reshape(vl,ndf*nel,1);
        
        ElemE = 1/2*vlres'*Me*vlres + 1/2*ulres'*Ke*ulres;
        ElemE;
        
end