% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       from mm10_form.f
% c
% c           Calculate the slip increment along system i  
% c           Dyson model
function slipinc = ...
    mm10_slipinc_dys(props, np1, n, stress, tt, alpha)

        dt = np1.tinc;
        T = np1.temp;
        k = props.boltzman;
        b = props.burgers; %m, from paper about IN738LC
        Na =  6.0221412927e23;
        R = k*Na;%8.314462175; % check units
%         time = np1.u(1);
        Qjd = props.c5; %kJ/mol, page 5 bottom left
        Ds0 = props.c6; %m^2/s, page 5 bottom left
        
        rho = tt(1);
        ms = np1.ms(1:6,alpha);
        rs = stress*ms; % tau^a
        t_m = abs(rs);
%         r_i = 12.5e-9; %nm, radius for gamma' in Nim90, page 5 middle left
        rhoT = rho; % maybe, not sure, based on page 463 of Power Plant
        G = props.mu_0; % MPA, from http://www.springworksutah.com/index.php/materials/nickel-base-alloy/niminoc-90.html
        a = 0.5; % maybe, not sure, based on page 463 of Power Plant
        %.5*G*b*sqrt(1e10)=1.03125e+06, which is "order unity",
        %like on page 4 middle left
        
%         Ts = 1000+273; % solvus temperature, Fig4 page 5, for Nim90, and page 4 middle right
%         C0 = 17*exp(-7250/Ts);
%         Ce = 17*exp(-7250/T); %(19), page 4 lower right
%         phi_p = (C0-Ce)/(0.23-Ce); %(18), page 4 middle right
%         r_p = (r_i^3 + exp(-3e5/R/T)*time)^(1/3); %(17), page 4 right column
        phi_p = np1.u(2);
        r_p = np1.u(3);
        % check units on time or the 3e5
        
        t_net = a*G*b*sqrt(rhoT); %
        np1.u(4) = t_net;
        t_D = max(t_m - t_net,0); %page 2, left column
        l_p = 1.6*r_p*(sqrt(pi/(4*phi_p))-1); % page 
        cjDs = Ds0*exp(-Qjd/R/T); % Power plant paper, page 459
        slipinc = dt * 1.6*rho*phi_p * (sqrt(pi/(4*phi_p))-1) * cjDs * sinh(t_D*b^2*l_p/k/T) * sign(rs);
        np1.u(5) = slipinc;
      
end
