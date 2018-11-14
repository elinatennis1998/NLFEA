% Compute tau_pass from Roters model
function tpass = mm10_tpass_halite(props, np1, n, rho, alpha)

% Set some properties
G = props.cp_033;
b = props.cp_040;
% c1 = [props.cp_001*ones(6,1);props.cp_002*ones(6,1);props.cp_003*ones(12,1)];
if alpha <= 6 || props.s_type == 6
      c1 = cp_001;
elseif alpha <= 12
      c1 = cp_002;
else
      c1 = cp_003;
end

% Get parallel dislocations
[~,rhoP] = mm10_rhoFP_halite(props, np1, n, rho, alpha);

% Evaluate the equation
tpass = c1.*G*b*sqrt(rhoP);