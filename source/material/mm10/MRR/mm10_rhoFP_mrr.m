% Compute parallel and forest densities for Roters model
function [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, rho, alpha)

% Set some properties
[G,H] = mm10_mrr_GH(props);

rhoF = G(alpha,:)*reshape(rho(1:12),12,1);
rhoP = H(alpha,:)*reshape(rho(1:12),12,1);