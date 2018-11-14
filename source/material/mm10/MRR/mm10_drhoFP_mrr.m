% Compute parallel and forest densities for Roters model
function [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, rho, alpha, beta)

% Set some properties
[G,H] = mm10_mrr_GH(props);

drhoF = G(alpha,beta);
drhoP = H(alpha,beta);