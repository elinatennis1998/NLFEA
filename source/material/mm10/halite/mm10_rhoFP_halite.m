% Compute parallel and forest densities for Roters model
function [rhoF,rhoP] = mm10_rhoFP_halite(props, np1, n, rho, alpha)

% Set some properties
[G,H] = mm10_halite_GH(props);

rhoF = G(alpha,:)*reshape(rho(1:24),24,1);
rhoP = H(alpha,:)*reshape(rho(1:24),24,1);