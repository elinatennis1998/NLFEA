% Compute tau_cut from Roters model
function tcut = mm10_tcut_halite(props, np1, n, rho, alpha)

% Set some properties
if props.s_type == 6
    Qslip = props.cp_025*ones(12,1);
    c1 = props.cp_001*ones(12,1);
    c2 = props.cp_004*ones(12,1);
elseif props.s_type == 11
    Qslip = [props.cp_025*ones(6,1);props.cp_026*ones(6,1);props.cp_027*ones(12,1)];
    c1 = [props.cp_001*ones(6,1);props.cp_002*ones(6,1);props.cp_003*ones(12,1)];
    c2 = [props.cp_004*ones(6,1);props.cp_005*ones(6,1);props.cp_006*ones(12,1)];
else
    error(' ');
end
b = props.burgers;
% Get forest dislocations
[rhoF,~] = mm10_rhoFP_halite(props, np1, n, rho, alpha);

% Evaluate the equation
tcut = Qslip(alpha)./(c1(alpha).*c2(alpha)*b^2).*sqrt(rhoF);