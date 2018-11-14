% Compute tau_pass from Roters model
function tpass = mm10_tpass_mrr(props, np1, n, rho, alpha)

% Set some properties
G = props.mu_0;
b = props.burgers;
c1 = props.c1;

% Get parallel dislocations
[~,rhoP] = mm10_rhoFP_mrr(props, np1, n, rho, alpha);

% Evaluate the equation
tpass = c1*G*b*sqrt(rhoP);