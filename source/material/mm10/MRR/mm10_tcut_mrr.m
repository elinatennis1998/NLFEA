% Compute tau_cut from Roters model
function tcut = mm10_tcut_mrr(props, np1, n, rho, alpha)

% Set some properties
Qslip = props.Qslip;
b = props.burgers;
c1 = props.c1;
c2 = props.c2;

% Get forest dislocations
[rhoF,~] = mm10_rhoFP_mrr(props, np1, n, rho, alpha);

% Evaluate the equation
tcut = Qslip/(c1*c2*b^2)*sqrt(rhoF);