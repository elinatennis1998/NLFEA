classdef crystal_state < handle % makes it act like pass-by-reference
% c     ****************************************************************
% c     *                                                              *
% c     *                 module mm10_defs                             *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 12/11/13                    *
% c     *                                                              *
% c     *    small module to hold crystal update data structs          *
% c     *                                                              *
% c     ****************************************************************
% c
%       module mm10_defs
%             implicit integer (a-z)
% $add param_def
% c           Massive list of properties
% max_slip_sys = 12;
% max_uhard = 20;
    properties
        R = zeros(3,3); Rp = zeros(3,3);
        stress = zeros(6,1); D = zeros(6,1);
        euler_angles = zeros(3,1);
        tau_l = zeros(maxparamsCP,1);
        slip_incs = zeros(maxparamsCP,1);
        gradFeinv = zeros(3,3,3);
        tangent = zeros(6,6);
        tau_tilde = zeros(maxparamsCP,1); temp = 0; tinc = 0;
        dg = 0; tau_v = 0; tau_y = 0;
        mu_harden = 0; work_inc = 0; p_work_inc = 0; p_strain_inc = 0;
        ms = zeros(6,maxparamsCP);
        qs = zeros(3,maxparamsCP); qc = zeros(3,maxparamsCP);
        u = zeros(maxparamsCP,1);
        step = 0; elem = 0; gp = 0;
        elaststrain = zeros(6,1);
        tt_rate = zeros(maxparamsCP,1);
        slip_incs_omar = zeros(12,1);     
        backstress_omar = zeros(12,1);
    end 
   
   methods
      function CS = crystal_state(step,elem) % Initialization method
          if nargin > 0
              CS.step  = step;
              CS.elem = elem;
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
end %classdef