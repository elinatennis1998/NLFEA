% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       Compute temperature dependent G modulus
function [G_110, G_100, G_111] = mm10_halite_Gtemp(T)
% Ran calculate shear modulus of halite
% Durand, M. A. (1936). "The temperature variation of the elastic moduli of
% NaCl, KCl and MgO." Physical Review 50(5): 449.
% T (K), p (Pa)
% calculate stiffness constant depending on temperature
c11_0 = 5.85e10;
c44_0 = 1.339e10;
A11 = 0.210;
A44 = 0.0601;
theta_80 = 320;
u = T/theta_80;
F11 = 0.2697 * u.^2 + 0.8293 * u - 0.1089;
F44 = 0.3676 * u.^2 + 0.7113 * u - 0.1092;
c11 = c11_0 * exp( -A11 * F11 );
c44 = c44_0 * exp( -A44 * F44 );
c12 = 1.17e10;

% shear modulus of {1 1 0}<1 -1 0>
G_110 = (c11 - c12) / 2;
% shear modulus of {1 0 0}<0  1 1>
G_100 = c44;
% shear modulus of {1 1 1}<1 -1 0>
G_111 = 3 * c44 .* (c11 - c12) ./ (c11 - c12 + 4 * c44);

end