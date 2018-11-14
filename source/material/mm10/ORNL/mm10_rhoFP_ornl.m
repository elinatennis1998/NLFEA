% Compute edge SSD density for slip system i fro AADP model
function [rhoF,rhoP] = mm10_rhoFP_ornl(props, np1, n, rho, alpha)

% Set some properties
[G,H] = mm10_ornl_GH(props);

rhoF = G(alpha,:)*reshape(rho(1:12),12,1);
rhoP = H(alpha,:)*reshape(rho(1:12),12,1);