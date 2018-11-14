% c
% c
% c     DYS:
% c
% c           Setup DYS hardening
function [props, np1, n] = mm10_setup_dys(props, np1, n)

% increment the total time
    time = n.u(1) + np1.tinc;
    np1.u(1) = time;
        T = np1.temp;
        k = props.boltzman;
        Na =  6.0221412927e23;
        R = k*Na;%8.314462175; % check units
%         time = np1.u(1);
        r_i = 12.5e-9; %nm, radius for gamma' in Nim90, page 5 middle left
        
        Ts = props.c3+273; % solvus temperature, Fig4 page 5, for Nim90, and page 4 middle right
        C0 = 17*exp(-7250/Ts);
        Ce = 17*exp(-7250/T); %(19), page 4 lower right
        phi_p = (C0-Ce)/(0.23-Ce); %(18), page 4 middle right
%         r_p = (r_i^3 + exp(-5.2e5/R/T)*time)^(1/3); %(17), page 4 right column
        r_p = (r_i^3 + exp(-props.c4/R/T)*time)^(1/3); %(17), page 4 right column
        % I changed to 5.2 to give more reasonable particle evolution
    np1.u(2) = phi_p;
    np1.u(3) = r_p;
% c     
%       return
end