classdef crystal_props < handle % makes it act like pass-by-reference
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
% This one lives in mm10
%       module mm10_defs
%             implicit integer (a-z)
% $add param_def
% c           Massive list of properties
% max_slip_sys = 12;
    properties
        rate_n= 0; tau_hat_y= 0; G_0_y= 0; burgers= 0; 
        p_v= 0; q_v= 0; boltzman= 0; theta_0= 0; eps_dot_0_v= 0; 
        eps_dot_0_y= 0; 
        p_y, q_y= 0;  
        tau_a= 0;  tau_hat_v= 0;  G_0_v= 0; 
        k_0= 0;  mu_0= 0;  D_0= 0;  T_0= 0;  tau_y= 0;  tau_v= 0;  voche_m= 0; 
        u1= 0;  u2= 0;  u3= 0;  u4= 0;  u5= 0;  u6= 0; 
        c1= 0; c2= 0; c3= 0; c4= 0; c5= 0; c6= 0; c7= 0; c8= 0;
        mfp= 0; v0=0; Qslip= 0; Qbulk = 0;
        g = zeros(3,3);
        ms = zeros(6,1);
        qs = zeros(3,1);
        ns = zeros(3,1);
        ni = zeros(maxparamsCP,3);
        bi = zeros(maxparamsCP,3);
        stiffness = zeros(6,6);
        init_angles = zeros(3,1);
        angle_type = 0;
        angle_convention = 0;
        nslip = 0;
        h_type = 0;
        s_type = 0;
        out = 0;
        num_hard = 0;
		usefsolve = 0;
        plugin = 0;
        miter = 0;
        cp_001 = 0; cp_002 = 0; cp_003 = 0; cp_004 = 0; cp_005 = 0;
        cp_006 = 0; cp_007 = 0; cp_008 = 0; cp_009 = 0; cp_010 = 0;
        cp_011 = 0; cp_012 = 0; cp_013 = 0; cp_014 = 0; cp_015 = 0;
        cp_016 = 0; cp_017 = 0; cp_018 = 0; cp_019 = 0; cp_020 = 0;
        cp_021 = 0; cp_022 = 0; cp_023 = 0; cp_024 = 0; cp_025 = 0;
        cp_026 = 0; cp_027 = 0; cp_028 = 0; cp_029 = 0; cp_030 = 0;
        cp_031 = 0; cp_032 = 0; cp_033 = 0; cp_034 = 0; cp_035 = 0;
        cp_036 = 0; cp_037 = 0; cp_038 = 0; cp_039 = 0; cp_040 = 0;
        cp_041 = 0; cp_042 = 0; cp_043 = 0; cp_044 = 0; cp_045 = 0;
        cp_046 = 0; cp_047 = 0; cp_048 = 0; cp_049 = 0; cp_050 = 0;
        cp_051 = 0; cp_052 = 0; cp_053 = 0; cp_054 = 0; cp_055 = 0;
        cp_056 = 0; cp_057 = 0; cp_058 = 0; cp_059 = 0; cp_060 = 0;
        cp_061 = 0; cp_062 = 0; cp_063 = 0; cp_064 = 0; cp_065 = 0;
        cp_066 = 0; cp_067 = 0; cp_068 = 0; cp_069 = 0; cp_070 = 0;
        cp_071 = 0; cp_072 = 0; cp_073 = 0; cp_074 = 0; cp_075 = 0;
        cp_076 = 0; cp_077 = 0; cp_078 = 0; cp_079 = 0; cp_080 = 0;
        cp_081 = 0; cp_082 = 0; cp_083 = 0; cp_084 = 0; cp_085 = 0;
        cp_086 = 0; cp_087 = 0; cp_088 = 0; cp_089 = 0; cp_090 = 0;
        cp_091 = 0; cp_092 = 0; cp_093 = 0; cp_094 = 0; cp_095 = 0;
        cp_096 = 0; cp_097 = 0; cp_098 = 0; cp_099 = 0; cp_100 = 0;
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