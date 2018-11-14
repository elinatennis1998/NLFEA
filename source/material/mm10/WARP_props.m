classdef WARP_props < handle % makes it act like pass-by-reference
% put other WARP_3D properties in here
% I need to look again at how Mark's properties flow into the code
    properties
        ncrystals = 0;
        angle_convention = 0;
        angle_type = 0;
        angles = zeros(3,1);
    end
    
     methods
        
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