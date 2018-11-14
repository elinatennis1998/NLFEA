classdef CP_all < handle % makes it act like pass-by-reference
    % c
    % c           Secondary type for CP crystal properties
    % max_slip_sys = 10;
    properties
        cp_prop = crystal;
        cp_stat = crystal_state;
        cp_other = WARP_props;
        cp_allTemps = 0; % default value for constant gp temperatures
    end %type
    
    methods
        function CP = CP_all(h_type) % Initialization method
            if nargin > 0
                CP.cp_prop = crystal(h_type);
            end
        end % Initialize
        
        % Make a copy of a handle object
        function new = copy(this)
            % Instantiate new object of the same class
            new = feval(class(this));
            
            % Copy all non-hidden properties
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
    end % methods
end % classdef