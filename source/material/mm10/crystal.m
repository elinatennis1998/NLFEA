classdef crystal < handle
    % From mod_crystals.f
    % This one comes from the input file
    properties
                  slip_type = 1;
% c                             1) fcc
% c                             2) bcc
% c                             3) single
                  elastic_type = 1;
% c                             1) isotropic
% c                             2) cubic
                  nslip = 0;
                  h_type = 1;
% c                             1) voche
% c                             2) mts
% c                             3) user
%                   e = 0; nu = 0; mu = 0; harden_n = 0; tau_a = 0;
%                                       tau_hat_y = 0; g_o_y = 0; b = 0; p_v = 0; q_v = 0;
%                                       p_y = 0; q_y = 0; boltz = 0; 
%                                       eps_dot_o_y = 0; t_o = 0;
%                                       theta_o = 0; tau_bar_o = 0;
%                                       tau_hat_v = 0; g_o_v = 0; theta_f = 0;
%                                       tau_t = 0; eps_dot_o_v = 0; k_o = 0;
%                                       mu_o = 0; D_o = 0; tau_y = 0; tau_v = 0;
%                                       voche_m = 0;
%                   u1 = 0; u2 = 0; u3 = 0; u4 = 0; u5 = 0; u6 = 0;
                  elast_stiff = zeros(6,6);
                  elast_flex = zeros(6,6);
                  ni = zeros(maxparamsCP,3);bi = zeros(maxparamsCP,3);
                  e = 30000;
                  nu = 0.3;
                  mu = 11538.5;
                  harden_n = 20;
                  tau_a = 0;
                  tau_hat_y = -1.0;
                  g_o_y = -1.0;
                  tau_hat_v = -1.0;
                  g_o_v = -1.0;
                  b = 3.5e-10;
                  p_v = 0.5;
                  q_v = 2;
                  p_y = 0.5;
                  q_y = 2;
                  boltz = 1.3806E-29;
                  eps_dot_o_y = 1.0E10;
                  eps_dot_o_v = 1.0E10;
                  t_o = 273.15;
                  theta_o = 57.7;
                  k_o = 0.0;
                  mu_o = 11538.5;
                  D_o = 0.0;
                  tau_y = 0.0;
                  tau_v = 0.0;
                  voche_m = 0.0;
                  R_omar = 0.0;
                  u1 = 0.0;
                  u2 = 0.0;
                  u3 = 0.0;
                  u4 = 0.0;
                  u5 = 0.0;
                  u6 = 0.0;
        c1= 0; c2= 0; c3= 0; c4= 0; c5= 0; c6= 0; c7= 0; c8= 0;
        mfp= 0; v0=0; Qslip= 0; Qbulk = 0;
                  usefsolve = 0;
                  plugin = 0;
                  miter = 50;
% c
                  valid = 0;
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
      function CP = crystal(h_type) % Initialization method
          if nargin > 0
              CP.h_type = h_type;
          end
      end % Initialize
      function c_array = finalize_new_crystal(c_array) % Initialization method
          if nargin > 0
              num = 1;
                  if (c_array(num).slip_type == 1) %then
                        c_array(num).nslip = 12;
                        f = 1/sqrt(2.0D0);
                        
                        c_array(num).bi(1,1)=0;
                        c_array(num).bi(1,2)=-f;
                        c_array(num).bi(1,3)=f;
                        c_array(num).bi(2,1)=f;
                        c_array(num).bi(2,2)=0;
                        c_array(num).bi(2,3)=-f;
                        c_array(num).bi(3,1)=-f;
                        c_array(num).bi(3,2)=f;
                        c_array(num).bi(3,3)=0;
                        c_array(num).bi(4,1)=0;
                        c_array(num).bi(4,2)=f;
                        c_array(num).bi(4,3)=f;
                        c_array(num).bi(5,1)=f;
                        c_array(num).bi(5,2)=0;
                        c_array(num).bi(5,3)=f;
                        c_array(num).bi(6,1)=f;
                        c_array(num).bi(6,2)=-f;
                        c_array(num).bi(6,3)=0;
                        c_array(num).bi(7,1)=0;
                        c_array(num).bi(7,2)=-f;
                        c_array(num).bi(7,3)=f;
                        c_array(num).bi(8,1)=-f;
                        c_array(num).bi(8,2)=0;
                        c_array(num).bi(8,3)=-f;
                        c_array(num).bi(9,1)=f;
                        c_array(num).bi(9,2)=f;
                        c_array(num).bi(9,3)=0;
                        c_array(num).bi(10,1)=0;
                        c_array(num).bi(10,2)=f;
                        c_array(num).bi(10,3)=f;
                        c_array(num).bi(11,1)=f;
                        c_array(num).bi(11,2)=0;
                        c_array(num).bi(11,3)=-f;
                        c_array(num).bi(12,1)=-f;
                        c_array(num).bi(12,2)=-f;
                        c_array(num).bi(12,3)=0;

                        f = 1/sqrt(3.0D0);
                        c_array(num).ni(1,1)=f;
                        c_array(num).ni(1,2)=f;
                        c_array(num).ni(1,3)=f;
                        c_array(num).ni(2,1)=f;
                        c_array(num).ni(2,2)=f;
                        c_array(num).ni(2,3)=f;
                        c_array(num).ni(3,1)=f;
                        c_array(num).ni(3,2)=f;
                        c_array(num).ni(3,3)=f;
                        c_array(num).ni(4,1)=-f;
                        c_array(num).ni(4,2)=-f;
                        c_array(num).ni(4,3)=f;
                        c_array(num).ni(5,1)=-f;
                        c_array(num).ni(5,2)=-f;
                        c_array(num).ni(5,3)=f;
                        c_array(num).ni(6,1)=-f;
                        c_array(num).ni(6,2)=-f;
                        c_array(num).ni(6,3)=f;
                        c_array(num).ni(7,1)=-f;
                        c_array(num).ni(7,2)=f;
                        c_array(num).ni(7,3)=f;
                        c_array(num).ni(8,1)=-f;
                        c_array(num).ni(8,2)=f;
                        c_array(num).ni(8,3)=f;
                        c_array(num).ni(9,1)=-f;
                        c_array(num).ni(9,2)=f;
                        c_array(num).ni(9,3)=f;
                        c_array(num).ni(10,1)=f;
                        c_array(num).ni(10,2)=-f;
                        c_array(num).ni(10,3)=f;
                        c_array(num).ni(11,1)=f;
                        c_array(num).ni(11,2)=-f;
                        c_array(num).ni(11,3)=f;
                        c_array(num).ni(12,1)=f;
                        c_array(num).ni(12,2)=-f;
                        c_array(num).ni(12,3)=f;
                  elseif (c_array(num).slip_type == 2) %then
                        c_array(num).nslip = 12;
                        f = 1/sqrt(2.0D0);
                        
                        c_array(num).ni(1,1)=0;
                        c_array(num).ni(1,2)=-f;
                        c_array(num).ni(1,3)=f;
                        c_array(num).ni(2,1)=f;
                        c_array(num).ni(2,2)=0;
                        c_array(num).ni(2,3)=-f;
                        c_array(num).ni(3,1)=-f;
                        c_array(num).ni(3,2)=f;
                        c_array(num).ni(3,3)=0;
                        c_array(num).ni(4,1)=0;
                        c_array(num).ni(4,2)=f;
                        c_array(num).ni(4,3)=f;
                        c_array(num).ni(5,1)=f;
                        c_array(num).ni(5,2)=0;
                        c_array(num).ni(5,3)=f;
                        c_array(num).ni(6,1)=f;
                        c_array(num).ni(6,2)=-f;
                        c_array(num).ni(6,3)=0;
                        c_array(num).ni(7,1)=0;
                        c_array(num).ni(7,2)=-f;
                        c_array(num).ni(7,3)=f;
                        c_array(num).ni(8,1)=-f;
                        c_array(num).ni(8,2)=0;
                        c_array(num).ni(8,3)=-f;
                        c_array(num).ni(9,1)=f;
                        c_array(num).ni(9,2)=f;
                        c_array(num).ni(9,3)=0;
                        c_array(num).ni(10,1)=0;
                        c_array(num).ni(10,2)=f;
                        c_array(num).ni(10,3)=f;
                        c_array(num).ni(11,1)=f;
                        c_array(num).ni(11,2)=0;
                        c_array(num).ni(11,3)=-f;
                        c_array(num).ni(12,1)=-f;
                        c_array(num).ni(12,2)=-f;
                        c_array(num).ni(12,3)=0;

                        f = 1/sqrt(3.0D0);
                        c_array(num).bi(1,1)=f;
                        c_array(num).bi(1,2)=f;
                        c_array(num).bi(1,3)=f;
                        c_array(num).bi(2,1)=f;
                        c_array(num).bi(2,2)=f;
                        c_array(num).bi(2,3)=f;
                        c_array(num).bi(3,1)=f;
                        c_array(num).bi(3,2)=f;
                        c_array(num).bi(3,3)=f;
                        c_array(num).bi(4,1)=-f;
                        c_array(num).bi(4,2)=-f;
                        c_array(num).bi(4,3)=f;
                        c_array(num).bi(5,1)=-f;
                        c_array(num).bi(5,2)=-f;
                        c_array(num).bi(5,3)=f;
                        c_array(num).bi(6,1)=-f;
                        c_array(num).bi(6,2)=-f;
                        c_array(num).bi(6,3)=f;
                        c_array(num).bi(7,1)=-f;
                        c_array(num).bi(7,2)=f;
                        c_array(num).bi(7,3)=f;
                        c_array(num).bi(8,1)=-f;
                        c_array(num).bi(8,2)=f;
                        c_array(num).bi(8,3)=f;
                        c_array(num).bi(9,1)=-f;
                        c_array(num).bi(9,2)=f;
                        c_array(num).bi(9,3)=f;
                        c_array(num).bi(10,1)=f;
                        c_array(num).bi(10,2)=-f;
                        c_array(num).bi(10,3)=f;
                        c_array(num).bi(11,1)=f;
                        c_array(num).bi(11,2)=-f;
                        c_array(num).bi(11,3)=f;
                        c_array(num).bi(12,1)=f;
                        c_array(num).bi(12,2)=-f;
                        c_array(num).bi(12,3)=f;
                  elseif (c_array(num).slip_type == 3) %then
                        c_array(num).nslip = 1;

                        c_array(num).bi(1,1) = 1.0;
                        c_array(num).bi(1,2) = 0.0;
                        c_array(num).bi(1,3) = 0.0;

                        c_array(num).ni(1,1) = 0.0;
                        c_array(num).ni(1,2) = 1.0;
                        c_array(num).ni(1,3) = 0.0;
                  elseif (c_array(num).slip_type == 4) % beam bending
                        c_array(num).nslip = 3;
                        
                        c_array(num).bi(1,1)=sqrt(3)/2;
                        c_array(num).bi(1,2)=.5;
                        c_array(num).bi(1,3)=0;
                        c_array(num).bi(2,1)=-sqrt(3)/2;
                        c_array(num).bi(2,2)=.5;
                        c_array(num).bi(2,3)=0;
                        c_array(num).bi(3,1)=0;
                        c_array(num).bi(3,2)=1;
                        c_array(num).bi(3,3)=0;

                        c_array(num).ni(1,1)=-.5;
                        c_array(num).ni(1,2)=sqrt(3)/2;
                        c_array(num).ni(1,3)=0;
                        c_array(num).ni(2,1)=.5;
                        c_array(num).ni(2,2)=sqrt(3)/2;
                        c_array(num).ni(2,3)=0;
                        c_array(num).ni(3,1)=-1;
                        c_array(num).ni(3,2)=0;
                        c_array(num).ni(3,3)=0;
                  elseif (c_array(num).slip_type == 5) % beam bending2
                        c_array(num).nslip = 1;
                        
                        c_array(num).bi(1,1)=sqrt(3)/2;
                        c_array(num).bi(1,2)=.5;
                        c_array(num).bi(1,3)=0;

                        c_array(num).ni(1,1)=-.5;
                        c_array(num).ni(1,2)=sqrt(3)/2;
                        c_array(num).ni(1,3)=0;
                  elseif (c_array(num).slip_type == 6) %then
                        c_array(num).nslip = 12;
                        f = 1/sqrt(2.0D0);
                        %edge
                        c_array(num).bi(1,1)=f;
                        c_array(num).bi(1,2)=-f;
                        c_array(num).bi(1,3)=0;
                        
                        c_array(num).bi(2,1)=f;
                        c_array(num).bi(2,2)=0;
                        c_array(num).bi(2,3)=-f;
                        
                        c_array(num).bi(3,1)=0;
                        c_array(num).bi(3,2)=f;
                        c_array(num).bi(3,3)=-f;
                        
                        c_array(num).bi(4,1)=f;
                        c_array(num).bi(4,2)=f;
                        c_array(num).bi(4,3)=0;
                        
                        c_array(num).bi(5,1)=f;
                        c_array(num).bi(5,2)=0;
                        c_array(num).bi(5,3)=f;
                        
                        c_array(num).bi(6,1)=0;
                        c_array(num).bi(6,2)=f;
                        c_array(num).bi(6,3)=-f;
                        
                        c_array(num).bi(7,1)=f;
                        c_array(num).bi(7,2)=f;
                        c_array(num).bi(7,3)=0;
                        
                        c_array(num).bi(8,1)=f;
                        c_array(num).bi(8,2)=0;
                        c_array(num).bi(8,3)=-f;
                        
                        c_array(num).bi(9,1)=0;
                        c_array(num).bi(9,2)=f;
                        c_array(num).bi(9,3)=f;
                        
                        c_array(num).bi(10,1)=f;
                        c_array(num).bi(10,2)=-f;
                        c_array(num).bi(10,3)=0;
                        
                        c_array(num).bi(11,1)=f;
                        c_array(num).bi(11,2)=0;
                        c_array(num).bi(11,3)=f;
                        
                        c_array(num).bi(12,1)=0;
                        c_array(num).bi(12,2)=f;
                        c_array(num).bi(12,3)=f;
%                         %screw
%                         f = 1/sqrt(6.0D0);
%                         c_array(num).bi(13,1)=f;
%                         c_array(num).bi(13,2)=f;
%                         c_array(num).bi(13,3)=-2.d0*f;
%                         
%                         c_array(num).bi(14,1)=f;
%                         c_array(num).bi(14,2)=-2.d0*f;
%                         c_array(num).bi(14,3)=f;
%                         
%                         c_array(num).bi(15,1)=-2.d0*f;
%                         c_array(num).bi(15,2)=f;
%                         c_array(num).bi(15,3)=f;
%                         
%                         c_array(num).bi(16,1)=f;
%                         c_array(num).bi(16,2)=-f;
%                         c_array(num).bi(16,3)=2.d0*f;
%                         
%                         c_array(num).bi(17,1)=f;
%                         c_array(num).bi(17,2)=2.d0*f;
%                         c_array(num).bi(17,3)=-f;
%                         
%                         c_array(num).bi(18,1)=2.d0*f;
%                         c_array(num).bi(18,2)=-f;
%                         c_array(num).bi(18,3)=f;
%                         
%                         f = 1/sqrt(3.0D0);
%                         c_array(num).bi(19,1)=f;
%                         c_array(num).bi(19,2)=f;
%                         c_array(num).bi(19,3)=f;
%                         
%                         c_array(num).bi(20,1)=-f;
%                         c_array(num).bi(20,2)=-f;
%                         c_array(num).bi(20,3)=-f;
%                         
%                         c_array(num).bi(21,1)=f;
%                         c_array(num).bi(21,2)=f;
%                         c_array(num).bi(21,3)=f;
%                         
%                         c_array(num).bi(22,1)=f;
%                         c_array(num).bi(22,2)=-f;
%                         c_array(num).bi(22,3)=-f;
%                         
%                         c_array(num).bi(23,1)=-f;
%                         c_array(num).bi(23,2)=f;
%                         c_array(num).bi(23,3)=f;
%                         
%                         c_array(num).bi(24,1)=f;
%                         c_array(num).bi(24,2)=f;
%                         c_array(num).bi(24,3)=-f;

                        
                        f = 1/sqrt(3.0D0);
                        %edge
                        c_array(num).ni(1,1)=f;
                        c_array(num).ni(1,2)=f;
                        c_array(num).ni(1,3)=f;
                        
                        c_array(num).ni(2,1)=-f;
                        c_array(num).ni(2,2)=-f;
                        c_array(num).ni(2,3)=-f;
                        
                        c_array(num).ni(3,1)=f;
                        c_array(num).ni(3,2)=f;
                        c_array(num).ni(3,3)=f;
                        
                        c_array(num).ni(4,1)=f;
                        c_array(num).ni(4,2)=-f;
                        c_array(num).ni(4,3)=-f;
                        
                        c_array(num).ni(5,1)=-f;
                        c_array(num).ni(5,2)=f;
                        c_array(num).ni(5,3)=f;
                        
                        c_array(num).ni(6,1)=f;
                        c_array(num).ni(6,2)=-f;
                        c_array(num).ni(6,3)=-f;
                        
                        c_array(num).ni(7,1)=f;
                        c_array(num).ni(7,2)=-f;
                        c_array(num).ni(7,3)=f;
                        
                        c_array(num).ni(8,1)=f;
                        c_array(num).ni(8,2)=-f;
                        c_array(num).ni(8,3)=f;
                        
                        c_array(num).ni(9,1)=-f;
                        c_array(num).ni(9,2)=f;
                        c_array(num).ni(9,3)=-f;
                        
                        c_array(num).ni(10,1)=-f;
                        c_array(num).ni(10,2)=-f;
                        c_array(num).ni(10,3)=f;
                        
                        c_array(num).ni(11,1)=-f;
                        c_array(num).ni(11,2)=-f;
                        c_array(num).ni(11,3)=f;
                        
                        c_array(num).ni(12,1)=f;
                        c_array(num).ni(12,2)=f;
                        c_array(num).ni(12,3)=-f;
%                         %screw
%                         c_array(num).ni(13,1)=-f;
%                         c_array(num).ni(13,2)=-f;
%                         c_array(num).ni(13,3)=-f;
%                         
%                         c_array(num).ni(14,1)=f;
%                         c_array(num).ni(14,2)=f;
%                         c_array(num).ni(14,3)=f;
%                         
%                         c_array(num).ni(15,1)=-f;
%                         c_array(num).ni(15,2)=-f;
%                         c_array(num).ni(15,3)=-f;
%                         
%                         c_array(num).ni(16,1)=-f;
%                         c_array(num).ni(16,2)=f;
%                         c_array(num).ni(16,3)=f;
%                         
%                         c_array(num).ni(17,1)=f;
%                         c_array(num).ni(17,2)=-f;
%                         c_array(num).ni(17,3)=-f;
%                         
%                         c_array(num).ni(18,1)=-f;
%                         c_array(num).ni(18,2)=-f;
%                         c_array(num).ni(18,3)=f;
%                         
%                         f = 1/sqrt(6.0D0);
%                         c_array(num).ni(19,1)=f;
%                         c_array(num).ni(19,2)=f;
%                         c_array(num).ni(19,3)=-2.d0*f;
%                         
%                         c_array(num).ni(20,1)=f;
%                         c_array(num).ni(20,2)=-2.d0*f;
%                         c_array(num).ni(20,3)=f;
%                         
%                         c_array(num).ni(21,1)=-2.d0*f;
%                         c_array(num).ni(21,2)=f;
%                         c_array(num).ni(21,3)=f;
%                         
%                         c_array(num).ni(22,1)=f;
%                         c_array(num).ni(22,2)=-f;
%                         c_array(num).ni(22,3)=2.d0*f;
%                         
%                         c_array(num).ni(23,1)=f;
%                         c_array(num).ni(23,2)=2.d0*f;
%                         c_array(num).ni(23,3)=-f;
%                         
%                         c_array(num).ni(24,1)=2.d0*f;
%                         c_array(num).ni(24,2)=-f;
%                         c_array(num).ni(24,3)=f;

                  elseif (c_array(num).slip_type == 7) %then % BCC 48 sys
                        c_array(num).nslip = 12;
                        z0 = 0.d0;
                        f2 = 1.d0/sqrt(2.0D0);
                        f3 = 1.d0/sqrt(3.0D0);
                        %{110}<111>
                        c_array(num).bi( 1,1)=-f3;
                        c_array(num).bi( 1,2)= f3;
                        c_array(num).bi( 1,3)= f3;
                        
                        c_array(num).bi( 2,1)= f3;
                        c_array(num).bi( 2,2)=-f3;
                        c_array(num).bi( 2,3)= f3;
                        
                        c_array(num).bi( 3,1)= f3;
                        c_array(num).bi( 3,2)= f3;
                        c_array(num).bi( 3,3)= f3;
                        
                        c_array(num).bi( 4,1)= f3;
                        c_array(num).bi( 4,2)= f3;
                        c_array(num).bi( 4,3)=-f3;
                        
                        c_array(num).bi( 5,1)= f3;
                        c_array(num).bi( 5,2)= f3;
                        c_array(num).bi( 5,3)=-f3;
                        
                        c_array(num).bi( 6,1)=-f3;
                        c_array(num).bi( 6,2)= f3;
                        c_array(num).bi( 6,3)= f3;
                        
                        c_array(num).bi( 7,1)= f3;
                        c_array(num).bi( 7,2)= f3;
                        c_array(num).bi( 7,3)= f3;
                        
                        c_array(num).bi( 8,1)= f3;
                        c_array(num).bi( 8,2)=-f3;
                        c_array(num).bi( 8,3)= f3;
                        
                        c_array(num).bi( 9,1)= f3;
                        c_array(num).bi( 9,2)= f3;
                        c_array(num).bi( 9,3)=-f3;
                        
                        c_array(num).bi(10,1)= f3;
                        c_array(num).bi(10,2)=-f3;
                        c_array(num).bi(10,3)= f3;
                        
                        c_array(num).bi(11,1)= f3;
                        c_array(num).bi(11,2)= f3;
                        c_array(num).bi(11,3)= f3;
                        
                        c_array(num).bi(12,1)=-f3;
                        c_array(num).bi(12,2)= f3;
                        c_array(num).bi(12,3)= f3;
                        
                        %{110}<111>
                        c_array(num).ni( 1,1)= f2;
                        c_array(num).ni( 1,2)= f2;
                        c_array(num).ni( 1,3)= z0;
                        
                        c_array(num).ni( 2,1)= f2;
                        c_array(num).ni( 2,2)= f2;
                        c_array(num).ni( 2,3)= z0;
                        
                        c_array(num).ni( 3,1)= f2;
                        c_array(num).ni( 3,2)=-f2;
                        c_array(num).ni( 3,3)= z0;
                        
                        c_array(num).ni( 4,1)= f2;
                        c_array(num).ni( 4,2)=-f2;
                        c_array(num).ni( 4,3)= z0;
                        
                        c_array(num).ni( 5,1)= f2;
                        c_array(num).ni( 5,2)= z0;
                        c_array(num).ni( 5,3)= f2;
                        
                        c_array(num).ni( 6,1)= f2;
                        c_array(num).ni( 6,2)= z0;
                        c_array(num).ni( 6,3)= f2;
                        
                        c_array(num).ni( 7,1)= f2;
                        c_array(num).ni( 7,2)= z0;
                        c_array(num).ni( 7,3)=-f2;
                        
                        c_array(num).ni( 8,1)= f2;
                        c_array(num).ni( 8,2)= z0;
                        c_array(num).ni( 8,3)=-f2;
                        
                        c_array(num).ni( 9,1)= z0;
                        c_array(num).ni( 9,2)= f2;
                        c_array(num).ni( 9,3)= f2;
                        
                        c_array(num).ni(10,1)= z0;
                        c_array(num).ni(10,2)= f2;
                        c_array(num).ni(10,3)= f2;
                        
                        c_array(num).ni(11,1)= z0;
                        c_array(num).ni(11,2)= f2;
                        c_array(num).ni(11,3)=-f2;
                        
                        c_array(num).ni(12,1)= z0;
                        c_array(num).ni(12,2)= f2;
                        c_array(num).ni(12,3)=-f2;

                  elseif (c_array(num).slip_type == 8) %then % BCC 48 sys
                        c_array(num).nslip = 48;
                        z0 = 0.d0;
                        f2 = 1.d0/sqrt(2.0D0);
                        f3 = 1.d0/sqrt(3.0D0);
                        f112 = 1.d0/sqrt(6.0D0);
                        f211 = 2.d0/sqrt(6.0D0);
                        f123 = 1.d0/sqrt(14.0D0);
                        f213 = 2.d0/sqrt(14.0D0);
                        f312 = 3.d0/sqrt(14.0D0);
                        %{110}<111>
                        c_array(num).bi( 1,1)=-f3;
                        c_array(num).bi( 1,2)= f3;
                        c_array(num).bi( 1,3)= f3;
                        
                        c_array(num).bi( 2,1)= f3;
                        c_array(num).bi( 2,2)=-f3;
                        c_array(num).bi( 2,3)= f3;
                        
                        c_array(num).bi( 3,1)= f3;
                        c_array(num).bi( 3,2)= f3;
                        c_array(num).bi( 3,3)= f3;
                        
                        c_array(num).bi( 4,1)= f3;
                        c_array(num).bi( 4,2)= f3;
                        c_array(num).bi( 4,3)=-f3;
                        
                        c_array(num).bi( 5,1)= f3;
                        c_array(num).bi( 5,2)= f3;
                        c_array(num).bi( 5,3)=-f3;
                        
                        c_array(num).bi( 6,1)=-f3;
                        c_array(num).bi( 6,2)= f3;
                        c_array(num).bi( 6,3)= f3;
                        
                        c_array(num).bi( 7,1)= f3;
                        c_array(num).bi( 7,2)= f3;
                        c_array(num).bi( 7,3)= f3;
                        
                        c_array(num).bi( 8,1)= f3;
                        c_array(num).bi( 8,2)=-f3;
                        c_array(num).bi( 8,3)= f3;
                        
                        c_array(num).bi( 9,1)= f3;
                        c_array(num).bi( 9,2)= f3;
                        c_array(num).bi( 9,3)=-f3;
                        
                        c_array(num).bi(10,1)= f3;
                        c_array(num).bi(10,2)=-f3;
                        c_array(num).bi(10,3)= f3;
                        
                        c_array(num).bi(11,1)= f3;
                        c_array(num).bi(11,2)= f3;
                        c_array(num).bi(11,3)= f3;
                        
                        c_array(num).bi(12,1)=-f3;
                        c_array(num).bi(12,2)= f3;
                        c_array(num).bi(12,3)= f3;
                        %{112}<111>
                        c_array(num).bi(13,1)= f3;
                        c_array(num).bi(13,2)= f3;
                        c_array(num).bi(13,3)=-f3;
                        
                        c_array(num).bi(14,1)= f3;
                        c_array(num).bi(14,2)=-f3;
                        c_array(num).bi(14,3)= f3;
                        
                        c_array(num).bi(15,1)=-f3;
                        c_array(num).bi(15,2)= f3;
                        c_array(num).bi(15,3)= f3;
                        
                        c_array(num).bi(16,1)= f3;
                        c_array(num).bi(16,2)= f3;
                        c_array(num).bi(16,3)= f3;
                        
                        c_array(num).bi(17,1)= f3;
                        c_array(num).bi(17,2)=-f3;
                        c_array(num).bi(17,3)= f3;
                        
                        c_array(num).bi(18,1)= f3;
                        c_array(num).bi(18,2)= f3;
                        c_array(num).bi(18,3)=-f3;
                        
                        c_array(num).bi(19,1)= f3;
                        c_array(num).bi(19,2)= f3;
                        c_array(num).bi(19,3)= f3;
                        
                        c_array(num).bi(20,1)=-f3;
                        c_array(num).bi(20,2)= f3;
                        c_array(num).bi(20,3)= f3;
                        
                        c_array(num).bi(21,1)=-f3;
                        c_array(num).bi(21,2)= f3;
                        c_array(num).bi(21,3)= f3;
                        
                        c_array(num).bi(22,1)= f3;
                        c_array(num).bi(22,2)= f3;
                        c_array(num).bi(22,3)= f3;
                        
                        c_array(num).bi(23,1)= f3;
                        c_array(num).bi(23,2)= f3;
                        c_array(num).bi(23,3)=-f3;
                        
                        c_array(num).bi(24,1)= f3;
                        c_array(num).bi(24,2)=-f3;
                        c_array(num).bi(24,3)= f3;
                        %{123}<111>
                        c_array(num).bi(25,1)= f3;
                        c_array(num).bi(25,2)= f3;
                        c_array(num).bi(25,3)=-f3;
                        
                        c_array(num).bi(26,1)= f3;
                        c_array(num).bi(26,2)=-f3;
                        c_array(num).bi(26,3)= f3;
                        
                        c_array(num).bi(27,1)=-f3;
                        c_array(num).bi(27,2)= f3;
                        c_array(num).bi(27,3)= f3;
                        
                        c_array(num).bi(28,1)= f3;
                        c_array(num).bi(28,2)= f3;
                        c_array(num).bi(28,3)= f3;
                        
                        c_array(num).bi(29,1)= f3;
                        c_array(num).bi(29,2)=-f3;
                        c_array(num).bi(29,3)= f3;
                        
                        c_array(num).bi(30,1)= f3;
                        c_array(num).bi(30,2)= f3;
                        c_array(num).bi(30,3)=-f3;
                        
                        c_array(num).bi(31,1)= f3;
                        c_array(num).bi(31,2)= f3;
                        c_array(num).bi(31,3)= f3;
                        
                        c_array(num).bi(32,1)=-f3;
                        c_array(num).bi(32,2)= f3;
                        c_array(num).bi(32,3)= f3;
                        
                        c_array(num).bi(33,1)= f3;
                        c_array(num).bi(33,2)= f3;
                        c_array(num).bi(33,3)=-f3;
                        
                        c_array(num).bi(34,1)= f3;
                        c_array(num).bi(34,2)=-f3;
                        c_array(num).bi(34,3)= f3;
                        
                        c_array(num).bi(35,1)=-f3;
                        c_array(num).bi(35,2)= f3;
                        c_array(num).bi(35,3)= f3;
                        
                        c_array(num).bi(36,1)= f3;
                        c_array(num).bi(36,2)= f3;
                        c_array(num).bi(36,3)= f3;
                        
                        c_array(num).bi(37,1)= f3;
                        c_array(num).bi(37,2)=-f3;
                        c_array(num).bi(37,3)= f3;
                        
                        c_array(num).bi(38,1)= f3;
                        c_array(num).bi(38,2)= f3;
                        c_array(num).bi(38,3)=-f3;
                        
                        c_array(num).bi(39,1)= f3;
                        c_array(num).bi(39,2)= f3;
                        c_array(num).bi(39,3)= f3;
                        
                        c_array(num).bi(40,1)=-f3;
                        c_array(num).bi(40,2)= f3;
                        c_array(num).bi(40,3)= f3;
                        
                        c_array(num).bi(41,1)=-f3;
                        c_array(num).bi(41,2)= f3;
                        c_array(num).bi(41,3)= f3;
                        
                        c_array(num).bi(42,1)= f3;
                        c_array(num).bi(42,2)= f3;
                        c_array(num).bi(42,3)= f3;
                        
                        c_array(num).bi(43,1)= f3;
                        c_array(num).bi(43,2)= f3;
                        c_array(num).bi(43,3)=-f3;
                        
                        c_array(num).bi(44,1)= f3;
                        c_array(num).bi(44,2)=-f3;
                        c_array(num).bi(44,3)= f3;
                        
                        c_array(num).bi(45,1)=-f3;
                        c_array(num).bi(45,2)= f3;
                        c_array(num).bi(45,3)= f3;
                        
                        c_array(num).bi(46,1)= f3;
                        c_array(num).bi(46,2)= f3;
                        c_array(num).bi(46,3)= f3;
                        
                        c_array(num).bi(47,1)= f3;
                        c_array(num).bi(47,2)= f3;
                        c_array(num).bi(47,3)=-f3;
                        
                        c_array(num).bi(48,1)= f3;
                        c_array(num).bi(48,2)=-f3;
                        c_array(num).bi(48,3)= f3;

                        
                        %{110}<111>
                        c_array(num).ni( 1,1)= f2;
                        c_array(num).ni( 1,2)= f2;
                        c_array(num).ni( 1,3)= z0;
                        
                        c_array(num).ni( 2,1)= f2;
                        c_array(num).ni( 2,2)= f2;
                        c_array(num).ni( 2,3)= z0;
                        
                        c_array(num).ni( 3,1)= f2;
                        c_array(num).ni( 3,2)=-f2;
                        c_array(num).ni( 3,3)= z0;
                        
                        c_array(num).ni( 4,1)= f2;
                        c_array(num).ni( 4,2)=-f2;
                        c_array(num).ni( 4,3)= z0;
                        
                        c_array(num).ni( 5,1)= f2;
                        c_array(num).ni( 5,2)= z0;
                        c_array(num).ni( 5,3)= f2;
                        
                        c_array(num).ni( 6,1)= f2;
                        c_array(num).ni( 6,2)= z0;
                        c_array(num).ni( 6,3)= f2;
                        
                        c_array(num).ni( 7,1)= f2;
                        c_array(num).ni( 7,2)= z0;
                        c_array(num).ni( 7,3)=-f2;
                        
                        c_array(num).ni( 8,1)= f2;
                        c_array(num).ni( 8,2)= z0;
                        c_array(num).ni( 8,3)=-f2;
                        
                        c_array(num).ni( 9,1)= z0;
                        c_array(num).ni( 9,2)= f2;
                        c_array(num).ni( 9,3)= f2;
                        
                        c_array(num).ni(10,1)= z0;
                        c_array(num).ni(10,2)= f2;
                        c_array(num).ni(10,3)= f2;
                        
                        c_array(num).ni(11,1)= z0;
                        c_array(num).ni(11,2)= f2;
                        c_array(num).ni(11,3)=-f2;
                        
                        c_array(num).ni(12,1)= z0;
                        c_array(num).ni(12,2)= f2;
                        c_array(num).ni(12,3)=-f2;
                        %{112}<111>
                        c_array(num).ni(13,1)= f112;
                        c_array(num).ni(13,2)= f112;
                        c_array(num).ni(13,3)= f211;
                        
                        c_array(num).ni(14,1)=-f112;
                        c_array(num).ni(14,2)= f112;
                        c_array(num).ni(14,3)= f211;
                        
                        c_array(num).ni(15,1)= f112;
                        c_array(num).ni(15,2)=-f112;
                        c_array(num).ni(15,3)= f211;
                        
                        c_array(num).ni(16,1)= f112;
                        c_array(num).ni(16,2)= f112;
                        c_array(num).ni(16,3)=-f211;
                        
                        c_array(num).ni(17,1)= f112;
                        c_array(num).ni(17,2)= f211;
                        c_array(num).ni(17,3)= f112;
                        
                        c_array(num).ni(18,1)=-f112;
                        c_array(num).ni(18,2)= f211;
                        c_array(num).ni(18,3)= f112;
                        
                        c_array(num).ni(19,1)= f112;
                        c_array(num).ni(19,2)=-f211;
                        c_array(num).ni(19,3)= f112;
                        
                        c_array(num).ni(20,1)= f112;
                        c_array(num).ni(20,2)= f211;
                        c_array(num).ni(20,3)=-f112;
                        
                        c_array(num).ni(21,1)= f211;
                        c_array(num).ni(21,2)= f112;
                        c_array(num).ni(21,3)= f112;
                        
                        c_array(num).ni(22,1)=-f211;
                        c_array(num).ni(22,2)= f112;
                        c_array(num).ni(22,3)= f112;
                        
                        c_array(num).ni(23,1)= f211;
                        c_array(num).ni(23,2)=-f112;
                        c_array(num).ni(23,3)= f112;
                        
                        c_array(num).ni(24,1)= f211;
                        c_array(num).ni(24,2)= f112;
                        c_array(num).ni(24,3)=-f112;
                        %{123}<111>
                        c_array(num).ni(25,1)= f123;
                        c_array(num).ni(25,2)= f213;
                        c_array(num).ni(25,3)= f312;
                        
                        c_array(num).ni(26,1)=-f123;
                        c_array(num).ni(26,2)= f213;
                        c_array(num).ni(26,3)= f312;
                        
                        c_array(num).ni(27,1)= f123;
                        c_array(num).ni(27,2)=-f213;
                        c_array(num).ni(27,3)= f312;
                        
                        c_array(num).ni(28,1)= f123;
                        c_array(num).ni(28,2)= f213;
                        c_array(num).ni(28,3)=-f312;
                        
                        c_array(num).ni(29,1)= f123;
                        c_array(num).ni(29,2)= f312;
                        c_array(num).ni(29,3)= f213;
                        
                        c_array(num).ni(30,1)=-f123;
                        c_array(num).ni(30,2)= f312;
                        c_array(num).ni(30,3)= f213;
                        
                        c_array(num).ni(31,1)= f123;
                        c_array(num).ni(31,2)=-f312;
                        c_array(num).ni(31,3)= f213;
                        
                        c_array(num).ni(32,1)= f123;
                        c_array(num).ni(32,2)= f312;
                        c_array(num).ni(32,3)=-f213;
                        
                        c_array(num).ni(33,1)= f213;
                        c_array(num).ni(33,2)= f123;
                        c_array(num).ni(33,3)= f312;
                        
                        c_array(num).ni(34,1)=-f213;
                        c_array(num).ni(34,2)= f123;
                        c_array(num).ni(34,3)= f312;
                        
                        c_array(num).ni(35,1)= f213;
                        c_array(num).ni(35,2)=-f123;
                        c_array(num).ni(35,3)= f312;
                        
                        c_array(num).ni(36,1)= f213;
                        c_array(num).ni(36,2)= f123;
                        c_array(num).ni(36,3)=-f312;

                        c_array(num).ni(37,1)= f213;
                        c_array(num).ni(37,2)= f312;
                        c_array(num).ni(37,3)= f123;
                        
                        c_array(num).ni(38,1)=-f213;
                        c_array(num).ni(38,2)= f312;
                        c_array(num).ni(38,3)= f123;
                        
                        c_array(num).ni(39,1)= f213;
                        c_array(num).ni(39,2)=-f312;
                        c_array(num).ni(39,3)= f123;
                        
                        c_array(num).ni(40,1)= f213;
                        c_array(num).ni(40,2)= f312;
                        c_array(num).ni(40,3)=-f123;
                        
                        c_array(num).ni(41,1)= f312;
                        c_array(num).ni(41,2)= f123;
                        c_array(num).ni(41,3)= f213;
                        
                        c_array(num).ni(42,1)=-f312;
                        c_array(num).ni(42,2)= f123;
                        c_array(num).ni(42,3)= f213;
                        
                        c_array(num).ni(43,1)= f312;
                        c_array(num).ni(43,2)=-f123;
                        c_array(num).ni(43,3)= f213;
                        
                        c_array(num).ni(44,1)= f312;
                        c_array(num).ni(44,2)= f123;
                        c_array(num).ni(44,3)=-f213;
                        
                        c_array(num).ni(45,1)= f312;
                        c_array(num).ni(45,2)= f213;
                        c_array(num).ni(45,3)= f123;
                        
                        c_array(num).ni(46,1)=-f312;
                        c_array(num).ni(46,2)= f213;
                        c_array(num).ni(46,3)= f123;
                        
                        c_array(num).ni(47,1)= f312;
                        c_array(num).ni(47,2)=-f213;
                        c_array(num).ni(47,3)= f123;
                        
                        c_array(num).ni(48,1)= f312;
                        c_array(num).ni(48,2)= f213;
                        c_array(num).ni(48,3)=-f123;
                        
                  elseif (c_array(num).slip_type == 9) %Ran HCP
                        c_array(num).nslip = 6;
                        % material constant of hcp
                        a = 0.295; c = 0.468; 
                        ac1 = sqrt(c^2+a^2);
                        ac2 = sqrt(4*c^2+3*a^2);
                        
                        % Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        % B1
                        c_array(num).bi(1,1)=0.5;
                        c_array(num).bi(1,2)=-sqrt(3)/2;
                        c_array(num).bi(1,3)=0;
                        % B2
                        c_array(num).bi(2,1)=0.5;
                        c_array(num).bi(2,2)=sqrt(3)/2;
                        c_array(num).bi(2,3)=0;
                        % B3
                        c_array(num).bi(3,1)=-1;
                        c_array(num).bi(3,2)=0;
                        c_array(num).bi(3,3)=0;
                        
                        % Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        % P1
                        c_array(num).bi(4,1)=1;
                        c_array(num).bi(4,2)=0;
                        c_array(num).bi(4,3)=0;
                        % P2
                        c_array(num).bi(5,1)=0.5;
                        c_array(num).bi(5,2)=sqrt(3)/2;
                        c_array(num).bi(5,3)=0;
                        % P3
                        c_array(num).bi(6,1)=-0.5;
                        c_array(num).bi(6,2)=sqrt(3)/2;
                        c_array(num).bi(6,3)=0;
                        
%                         % Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
%                         % R1
%                         c_array(num).bi(7,1)=1;
%                         c_array(num).bi(7,2)=0;
%                         c_array(num).bi(7,3)=0;
%                         % R2
%                         c_array(num).bi(8,1)=0.5;
%                         c_array(num).bi(8,2)=sqrt(3)/2;
%                         c_array(num).bi(8,3)=0;
%                         % R3
%                         c_array(num).bi(9,1)=-0.5;
%                         c_array(num).bi(9,2)=sqrt(3)/2;
%                         c_array(num).bi(9,3)=0;
%                         % R4
%                         c_array(num).bi(10,1)=-1;
%                         c_array(num).bi(10,2)=0;
%                         c_array(num).bi(10,3)=0;
%                         % R5
%                         c_array(num).bi(11,1)=-0.5;
%                         c_array(num).bi(11,2)=-sqrt(3)/2;
%                         c_array(num).bi(11,3)=0;
%                         % R6
%                         c_array(num).bi(12,1)=0.5;
%                         c_array(num).bi(12,2)=-sqrt(3)/2;
%                         c_array(num).bi(12,3)=0;
%                         
%                         % Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
%                         % R7
%                         c_array(num).bi(13,1)=a/(2*ac1);
%                         c_array(num).bi(13,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(13,3)=c/ac1;
%                         % R8
%                         c_array(num).bi(14,1)=-a/(2*ac1);
%                         c_array(num).bi(14,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(14,3)=c/ac1;
%                         % R9
%                         c_array(num).bi(15,1)=-a/ac1;
%                         c_array(num).bi(15,2)=0;
%                         c_array(num).bi(15,3)=c/ac1;
%                         % R10
%                         c_array(num).bi(16,1)=-a/(2*ac1);
%                         c_array(num).bi(16,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(16,3)=c/ac1;
%                         % R11
%                         c_array(num).bi(17,1)=a/(2*ac1);
%                         c_array(num).bi(17,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(17,3)=c/ac1;
%                         % R12
%                         c_array(num).bi(18,1)=a/ac1;
%                         c_array(num).bi(18,2)=0;
%                         c_array(num).bi(18,3)=c/ac1;
%                         % R13
%                         c_array(num).bi(19,1)=-a/(2*ac1);
%                         c_array(num).bi(19,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(19,3)=c/ac1;
%                         % R14
%                         c_array(num).bi(20,1)=-a/ac1;
%                         c_array(num).bi(20,2)=0;
%                         c_array(num).bi(20,3)=c/ac1;
%                         % R15
%                         c_array(num).bi(21,1)=-a/(2*ac1);
%                         c_array(num).bi(21,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(21,3)=c/ac1;
%                         % R16
%                         c_array(num).bi(22,1)=a/(2*ac1);
%                         c_array(num).bi(22,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(22,3)=c/ac1;
%                         % R17
%                         c_array(num).bi(23,1)=a/ac1;
%                         c_array(num).bi(23,2)=0;
%                         c_array(num).bi(23,3)=c/ac1;
%                         % R18
%                         c_array(num).bi(24,1)=a/(2*ac1);
%                         c_array(num).bi(24,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(24,3)=c/ac1;
% 
%                         % Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
%                         % R19
%                         c_array(num).bi(25,1)=-a/(2*ac1);
%                         c_array(num).bi(25,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(25,3)=c/ac1;
%                         % R20
%                         c_array(num).bi(26,1)=-a/ac1;
%                         c_array(num).bi(26,2)=0;
%                         c_array(num).bi(26,3)=c/ac1;
%                         % R21
%                         c_array(num).bi(27,1)=-a/(2*ac1);
%                         c_array(num).bi(27,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(27,3)=c/ac1;
%                         % R22
%                         c_array(num).bi(28,1)=a/(2*ac1);
%                         c_array(num).bi(28,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(28,3)=c/ac1;
%                         % R23
%                         c_array(num).bi(29,1)=a/ac1;
%                         c_array(num).bi(29,2)=0;
%                         c_array(num).bi(29,3)=c/ac1;
%                         % R24
%                         c_array(num).bi(30,1)=a/(2*ac1);
%                         c_array(num).bi(30,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(30,3)=c/ac1;
% -------------------------------------------------------------------------
                        % Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        % B1
                        c_array(num).ni(1,1)=0;
                        c_array(num).ni(1,2)=0;
                        c_array(num).ni(1,3)=1;
                        % B2
                        c_array(num).ni(2,1)=0;
                        c_array(num).ni(2,2)=0;
                        c_array(num).ni(2,3)=1;
                        % B3
                        c_array(num).ni(3,1)=0;
                        c_array(num).ni(3,2)=0;
                        c_array(num).ni(3,3)=1;
                        
                        % Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        % P1
                        c_array(num).ni(4,1)=0;
                        c_array(num).ni(4,2)=1;
                        c_array(num).ni(4,3)=0;
                        % P2
                        c_array(num).ni(5,1)=-sqrt(3)/2;
                        c_array(num).ni(5,2)=0.5;
                        c_array(num).ni(5,3)=0;
                        % P3
                        c_array(num).ni(6,1)=-sqrt(3)/2;
                        c_array(num).ni(6,2)=-0.5;
                        c_array(num).ni(6,3)=0;
                        
%                         % Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
%                         % R1
%                         c_array(num).ni(7,1)=0;
%                         c_array(num).ni(7,2)=-2*c/ac2;
%                         c_array(num).ni(7,3)=sqrt(3)*a/ac2;
%                         % R2
%                         c_array(num).ni(8,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(8,2)=-c/ac2;
%                         c_array(num).ni(8,3)=sqrt(3)*a/ac2;
%                         % R3
%                         c_array(num).ni(9,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(9,2)=c/ac2;
%                         c_array(num).ni(9,3)=sqrt(3)*a/ac2;
%                         % R4
%                         c_array(num).ni(10,1)=0;
%                         c_array(num).ni(10,2)=2*c/ac2;
%                         c_array(num).ni(10,3)=sqrt(3)*a/ac2;
%                         % R5
%                         c_array(num).ni(11,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(11,2)=c/ac2;
%                         c_array(num).ni(11,3)=sqrt(3)*a/ac2;
%                         % R6
%                         c_array(num).ni(12,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(12,2)=-c/ac2;
%                         c_array(num).ni(12,3)=sqrt(3)*a/ac2;
%                         
%                         % Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
%                         % R7
%                         c_array(num).ni(13,1)=0;
%                         c_array(num).ni(13,2)=-2*c/ac2;
%                         c_array(num).ni(13,3)=sqrt(3)*a/ac2;
%                         % R8
%                         c_array(num).ni(14,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(14,2)=-c/ac2;
%                         c_array(num).ni(14,3)=sqrt(3)*a/ac2;
%                         % R9
%                         c_array(num).ni(15,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(15,2)=c/ac2;
%                         c_array(num).ni(15,3)=sqrt(3)*a/ac2;
%                         % R10
%                         c_array(num).ni(16,1)=0;
%                         c_array(num).ni(16,2)=2*c/ac2;
%                         c_array(num).ni(16,3)=sqrt(3)*a/ac2;
%                         % R11
%                         c_array(num).ni(17,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(17,2)=c/ac2;
%                         c_array(num).ni(17,3)=sqrt(3)*a/ac2;
%                         % R12
%                         c_array(num).ni(18,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(18,2)=-c/ac2;
%                         c_array(num).ni(18,3)=sqrt(3)*a/ac2;
%                         % R13
%                         c_array(num).ni(19,1)=0;
%                         c_array(num).ni(19,2)=-2*c/ac2;
%                         c_array(num).ni(19,3)=sqrt(3)*a/ac2;
%                         % R14
%                         c_array(num).ni(20,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(20,2)=-c/ac2;
%                         c_array(num).ni(20,3)=sqrt(3)*a/ac2;
%                         % R15
%                         c_array(num).ni(21,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(21,2)=c/ac2;
%                         c_array(num).ni(21,3)=sqrt(3)*a/ac2;
%                         % R16
%                         c_array(num).ni(22,1)=0;
%                         c_array(num).ni(22,2)=2*c/ac2;
%                         c_array(num).ni(22,3)=sqrt(3)*a/ac2;
%                         % R17
%                         c_array(num).ni(23,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(23,2)=c/ac2;
%                         c_array(num).ni(23,3)=sqrt(3)*a/ac2;
%                         % R18
%                         c_array(num).ni(24,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(24,2)=-c/ac2;
%                         c_array(num).ni(24,3)=sqrt(3)*a/ac2;
% 
%                         % Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
%                         % R19
%                         c_array(num).ni(25,1)=c/(2*ac1);
%                         c_array(num).ni(25,2)=-sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(25,3)=a/ac1;
%                         % R20
%                         c_array(num).ni(26,1)=c/ac1;
%                         c_array(num).ni(26,2)=0;
%                         c_array(num).ni(26,3)=a/ac1;
%                         % R21
%                         c_array(num).ni(27,1)=c/(2*ac1);
%                         c_array(num).ni(27,2)=sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(27,3)=a/ac1;
%                         % R22
%                         c_array(num).ni(28,1)=-c/(2*ac1);
%                         c_array(num).ni(28,2)=sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(28,3)=a/ac1;
%                         % R23
%                         c_array(num).ni(29,1)=-c/ac1;
%                         c_array(num).ni(29,2)=0;
%                         c_array(num).ni(29,3)=a/ac1;
%                         % R24
%                         c_array(num).ni(30,1)=-c/(2*ac1);
%                         c_array(num).ni(30,2)=-sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(30,3)=a/ac1;
                        
                  elseif (c_array(num).slip_type == 10) % HCP: basal, prismatic, 1st-order pyramidal
                        c_array(num).nslip = 6+12;
                        % material constant of hcp
                        a = 0.295; c = 0.468; 
                        ac1 = sqrt(c^2+a^2);
                        ac2 = sqrt(4*c^2+3*a^2);
                        
                        % Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        % B1
                        c_array(num).bi(1,1)=0.5;
                        c_array(num).bi(1,2)=-sqrt(3)/2;
                        c_array(num).bi(1,3)=0;
                        % B2
                        c_array(num).bi(2,1)=0.5;
                        c_array(num).bi(2,2)=sqrt(3)/2;
                        c_array(num).bi(2,3)=0;
                        % B3
                        c_array(num).bi(3,1)=-1;
                        c_array(num).bi(3,2)=0;
                        c_array(num).bi(3,3)=0;
                        
                        % Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        % P1
                        c_array(num).bi(4,1)=1;
                        c_array(num).bi(4,2)=0;
                        c_array(num).bi(4,3)=0;
                        % P2
                        c_array(num).bi(5,1)=0.5;
                        c_array(num).bi(5,2)=sqrt(3)/2;
                        c_array(num).bi(5,3)=0;
                        % P3
                        c_array(num).bi(6,1)=-0.5;
                        c_array(num).bi(6,2)=sqrt(3)/2;
                        c_array(num).bi(6,3)=0;
                        
%                         % Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
%                         % R1
%                         c_array(num).bi(7,1)=1;
%                         c_array(num).bi(7,2)=0;
%                         c_array(num).bi(7,3)=0;
%                         % R2
%                         c_array(num).bi(8,1)=0.5;
%                         c_array(num).bi(8,2)=sqrt(3)/2;
%                         c_array(num).bi(8,3)=0;
%                         % R3
%                         c_array(num).bi(9,1)=-0.5;
%                         c_array(num).bi(9,2)=sqrt(3)/2;
%                         c_array(num).bi(9,3)=0;
%                         % R4
%                         c_array(num).bi(10,1)=-1;
%                         c_array(num).bi(10,2)=0;
%                         c_array(num).bi(10,3)=0;
%                         % R5
%                         c_array(num).bi(11,1)=-0.5;
%                         c_array(num).bi(11,2)=-sqrt(3)/2;
%                         c_array(num).bi(11,3)=0;
%                         % R6
%                         c_array(num).bi(12,1)=0.5;
%                         c_array(num).bi(12,2)=-sqrt(3)/2;
%                         c_array(num).bi(12,3)=0;
%                         
                        % Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
                        % R7
                        c_array(num).bi(13-6,1)=a/(2*ac1);
                        c_array(num).bi(13-6,2)=sqrt(3)*a/(2*ac1);
                        c_array(num).bi(13-6,3)=c/ac1;
                        % R8
                        c_array(num).bi(14-6,1)=-a/(2*ac1);
                        c_array(num).bi(14-6,2)=sqrt(3)*a/(2*ac1);
                        c_array(num).bi(14-6,3)=c/ac1;
                        % R9
                        c_array(num).bi(15-6,1)=-a/ac1;
                        c_array(num).bi(15-6,2)=0;
                        c_array(num).bi(15-6,3)=c/ac1;
                        % R10
                        c_array(num).bi(16-6,1)=-a/(2*ac1);
                        c_array(num).bi(16-6,2)=-sqrt(3)*a/(2*ac1);
                        c_array(num).bi(16-6,3)=c/ac1;
                        % R11
                        c_array(num).bi(17-6,1)=a/(2*ac1);
                        c_array(num).bi(17-6,2)=-sqrt(3)*a/(2*ac1);
                        c_array(num).bi(17-6,3)=c/ac1;
                        % R12
                        c_array(num).bi(18-6,1)=a/ac1;
                        c_array(num).bi(18-6,2)=0;
                        c_array(num).bi(18-6,3)=c/ac1;
                        % R13
                        c_array(num).bi(19-6,1)=-a/(2*ac1);
                        c_array(num).bi(19-6,2)=sqrt(3)*a/(2*ac1);
                        c_array(num).bi(19-6,3)=c/ac1;
                        % R14
                        c_array(num).bi(20-6,1)=-a/ac1;
                        c_array(num).bi(20-6,2)=0;
                        c_array(num).bi(20-6,3)=c/ac1;
                        % R15
                        c_array(num).bi(21-6,1)=-a/(2*ac1);
                        c_array(num).bi(21-6,2)=-sqrt(3)*a/(2*ac1);
                        c_array(num).bi(21-6,3)=c/ac1;
                        % R16
                        c_array(num).bi(22-6,1)=a/(2*ac1);
                        c_array(num).bi(22-6,2)=-sqrt(3)*a/(2*ac1);
                        c_array(num).bi(22-6,3)=c/ac1;
                        % R17
                        c_array(num).bi(23-6,1)=a/ac1;
                        c_array(num).bi(23-6,2)=0;
                        c_array(num).bi(23-6,3)=c/ac1;
                        % R18
                        c_array(num).bi(24-6,1)=a/(2*ac1);
                        c_array(num).bi(24-6,2)=sqrt(3)*a/(2*ac1);
                        c_array(num).bi(24-6,3)=c/ac1;
% 
%                         % Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
%                         % R19
%                         c_array(num).bi(25,1)=-a/(2*ac1);
%                         c_array(num).bi(25,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(25,3)=c/ac1;
%                         % R20
%                         c_array(num).bi(26,1)=-a/ac1;
%                         c_array(num).bi(26,2)=0;
%                         c_array(num).bi(26,3)=c/ac1;
%                         % R21
%                         c_array(num).bi(27,1)=-a/(2*ac1);
%                         c_array(num).bi(27,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(27,3)=c/ac1;
%                         % R22
%                         c_array(num).bi(28,1)=a/(2*ac1);
%                         c_array(num).bi(28,2)=-sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(28,3)=c/ac1;
%                         % R23
%                         c_array(num).bi(29,1)=a/ac1;
%                         c_array(num).bi(29,2)=0;
%                         c_array(num).bi(29,3)=c/ac1;
%                         % R24
%                         c_array(num).bi(30,1)=a/(2*ac1);
%                         c_array(num).bi(30,2)=sqrt(3)*a/(2*ac1);
%                         c_array(num).bi(30,3)=c/ac1;
% -------------------------------------------------------------------------
                        % Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        % B1
                        c_array(num).ni(1,1)=0;
                        c_array(num).ni(1,2)=0;
                        c_array(num).ni(1,3)=1;
                        % B2
                        c_array(num).ni(2,1)=0;
                        c_array(num).ni(2,2)=0;
                        c_array(num).ni(2,3)=1;
                        % B3
                        c_array(num).ni(3,1)=0;
                        c_array(num).ni(3,2)=0;
                        c_array(num).ni(3,3)=1;
                        
                        % Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        % P1
                        c_array(num).ni(4,1)=0;
                        c_array(num).ni(4,2)=1;
                        c_array(num).ni(4,3)=0;
                        % P2
                        c_array(num).ni(5,1)=-sqrt(3)/2;
                        c_array(num).ni(5,2)=0.5;
                        c_array(num).ni(5,3)=0;
                        % P3
                        c_array(num).ni(6,1)=-sqrt(3)/2;
                        c_array(num).ni(6,2)=-0.5;
                        c_array(num).ni(6,3)=0;
                        
%                         % Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
%                         % R1
%                         c_array(num).ni(7,1)=0;
%                         c_array(num).ni(7,2)=-2*c/ac2;
%                         c_array(num).ni(7,3)=sqrt(3)*a/ac2;
%                         % R2
%                         c_array(num).ni(8,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(8,2)=-c/ac2;
%                         c_array(num).ni(8,3)=sqrt(3)*a/ac2;
%                         % R3
%                         c_array(num).ni(9,1)=sqrt(3)*c/ac2;
%                         c_array(num).ni(9,2)=c/ac2;
%                         c_array(num).ni(9,3)=sqrt(3)*a/ac2;
%                         % R4
%                         c_array(num).ni(10,1)=0;
%                         c_array(num).ni(10,2)=2*c/ac2;
%                         c_array(num).ni(10,3)=sqrt(3)*a/ac2;
%                         % R5
%                         c_array(num).ni(11,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(11,2)=c/ac2;
%                         c_array(num).ni(11,3)=sqrt(3)*a/ac2;
%                         % R6
%                         c_array(num).ni(12,1)=-sqrt(3)*c/ac2;
%                         c_array(num).ni(12,2)=-c/ac2;
%                         c_array(num).ni(12,3)=sqrt(3)*a/ac2;
%                         
                        % Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
                        % R7
                        c_array(num).ni(13-6,1)=0;
                        c_array(num).ni(13-6,2)=-2*c/ac2;
                        c_array(num).ni(13-6,3)=sqrt(3)*a/ac2;
                        % R8
                        c_array(num).ni(14-6,1)=sqrt(3)*c/ac2;
                        c_array(num).ni(14-6,2)=-c/ac2;
                        c_array(num).ni(14-6,3)=sqrt(3)*a/ac2;
                        % R9
                        c_array(num).ni(15-6,1)=sqrt(3)*c/ac2;
                        c_array(num).ni(15-6,2)=c/ac2;
                        c_array(num).ni(15-6,3)=sqrt(3)*a/ac2;
                        % R10
                        c_array(num).ni(16-6,1)=0;
                        c_array(num).ni(16-6,2)=2*c/ac2;
                        c_array(num).ni(16-6,3)=sqrt(3)*a/ac2;
                        % R11
                        c_array(num).ni(17-6,1)=-sqrt(3)*c/ac2;
                        c_array(num).ni(17-6,2)=c/ac2;
                        c_array(num).ni(17-6,3)=sqrt(3)*a/ac2;
                        % R12
                        c_array(num).ni(18-6,1)=-sqrt(3)*c/ac2;
                        c_array(num).ni(18-6,2)=-c/ac2;
                        c_array(num).ni(18-6,3)=sqrt(3)*a/ac2;
                        % R13
                        c_array(num).ni(19-6,1)=0;
                        c_array(num).ni(19-6,2)=-2*c/ac2;
                        c_array(num).ni(19-6,3)=sqrt(3)*a/ac2;
                        % R14
                        c_array(num).ni(20-6,1)=sqrt(3)*c/ac2;
                        c_array(num).ni(20-6,2)=-c/ac2;
                        c_array(num).ni(20-6,3)=sqrt(3)*a/ac2;
                        % R15
                        c_array(num).ni(21-6,1)=sqrt(3)*c/ac2;
                        c_array(num).ni(21-6,2)=c/ac2;
                        c_array(num).ni(21-6,3)=sqrt(3)*a/ac2;
                        % R16
                        c_array(num).ni(22-6,1)=0;
                        c_array(num).ni(22-6,2)=2*c/ac2;
                        c_array(num).ni(22-6,3)=sqrt(3)*a/ac2;
                        % R17
                        c_array(num).ni(23-6,1)=-sqrt(3)*c/ac2;
                        c_array(num).ni(23-6,2)=c/ac2;
                        c_array(num).ni(23-6,3)=sqrt(3)*a/ac2;
                        % R18
                        c_array(num).ni(24-6,1)=-sqrt(3)*c/ac2;
                        c_array(num).ni(24-6,2)=-c/ac2;
                        c_array(num).ni(24-6,3)=sqrt(3)*a/ac2;
% 
%                         % Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
%                         % R19
%                         c_array(num).ni(25,1)=c/(2*ac1);
%                         c_array(num).ni(25,2)=-sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(25,3)=a/ac1;
%                         % R20
%                         c_array(num).ni(26,1)=c/ac1;
%                         c_array(num).ni(26,2)=0;
%                         c_array(num).ni(26,3)=a/ac1;
%                         % R21
%                         c_array(num).ni(27,1)=c/(2*ac1);
%                         c_array(num).ni(27,2)=sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(27,3)=a/ac1;
%                         % R22
%                         c_array(num).ni(28,1)=-c/(2*ac1);
%                         c_array(num).ni(28,2)=sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(28,3)=a/ac1;
%                         % R23
%                         c_array(num).ni(29,1)=-c/ac1;
%                         c_array(num).ni(29,2)=0;
%                         c_array(num).ni(29,3)=a/ac1;
%                         % R24
%                         c_array(num).ni(30,1)=-c/(2*ac1);
%                         c_array(num).ni(30,2)=-sqrt(3)*c/(2*ac1);
%                         c_array(num).ni(30,3)=a/ac1;

%--------------------------------Ran Halite Start--------------------------
                  elseif (c_array(num).slip_type == 11) %then
                        % initialize material parameters for MRR model
                        if isequal(size(c_array(num).c1),[3,1])
                        HaliteMode = [repmat([1,0,0],6,1); 
                                      repmat([0,1,0],6,1); 
                                      repmat([0,0,1],12,1)];
                        c_array(num).c1 = HaliteMode*c_array(num).c1;
                        c_array(num).c2 = HaliteMode*c_array(num).c2;
                        c_array(num).c3 = HaliteMode*c_array(num).c3;
                        c_array(num).c4 = HaliteMode*c_array(num).c4;
                        c_array(num).c5 = HaliteMode*c_array(num).c5;
                        c_array(num).c6 = HaliteMode*c_array(num).c6;
                        c_array(num).c7 = HaliteMode*c_array(num).c7;
                        c_array(num).c8 = HaliteMode*c_array(num).c8;
                        c_array(num).Qslip = HaliteMode*c_array(num).Qslip;
                        c_array(num).Qbulk = HaliteMode*c_array(num).Qbulk;
                        end
                        c_array(num).nslip = 24;
                        g = 1/sqrt(2.0D0);
                        % slip direction
        %-------------------{1 1 0}<1 -1 0>------------------------
                        c_array(num).bi(1,1)=g;
                        c_array(num).bi(1,2)=-g;
                        c_array(num).bi(1,3)=0;
                        
                        c_array(num).bi(2,1)=g;
                        c_array(num).bi(2,2)=g;
                        c_array(num).bi(2,3)=0;
                        
                        c_array(num).bi(3,1)=g;
                        c_array(num).bi(3,2)=0;
                        c_array(num).bi(3,3)=-g;
                        
                        c_array(num).bi(4,1)=g;
                        c_array(num).bi(4,2)=0;
                        c_array(num).bi(4,3)=g;
                        
                        c_array(num).bi(5,1)=0;
                        c_array(num).bi(5,2)=g;
                        c_array(num).bi(5,3)=-g;
                        
                        c_array(num).bi(6,1)=0;
                        c_array(num).bi(6,2)=g;
                        c_array(num).bi(6,3)=g;
        %-------------------{1 0 0}<0  1 1>------------------------
                        c_array(num).bi(7,1)=0;
                        c_array(num).bi(7,2)=g;
                        c_array(num).bi(7,3)=g;
                        
                        c_array(num).bi(8,1)=0;
                        c_array(num).bi(8,2)=g;
                        c_array(num).bi(8,3)=-g;
                        
                        c_array(num).bi(9,1)=g;
                        c_array(num).bi(9,2)=0;
                        c_array(num).bi(9,3)=g;
                        
                        c_array(num).bi(10,1)=g;
                        c_array(num).bi(10,2)=0;
                        c_array(num).bi(10,3)=-g;
                        
                        c_array(num).bi(11,1)=g;
                        c_array(num).bi(11,2)=g;
                        c_array(num).bi(11,3)=0;
                        
                        c_array(num).bi(12,1)=g;
                        c_array(num).bi(12,2)=-g;
                        c_array(num).bi(12,3)=0;
        %-------------------{1 1 1}<1 -1 0>------------------------
                        c_array(num).bi(13,1)=g;
                        c_array(num).bi(13,2)=-g;
                        c_array(num).bi(13,3)=0;
                        
                        c_array(num).bi(14,1)=g;
                        c_array(num).bi(14,2)=0;
                        c_array(num).bi(14,3)=-g;
                        
                        c_array(num).bi(15,1)=0;
                        c_array(num).bi(15,2)=g;
                        c_array(num).bi(15,3)=-g;
                        
                        c_array(num).bi(16,1)=g;
                        c_array(num).bi(16,2)=g;
                        c_array(num).bi(16,3)=0;
                        
                        c_array(num).bi(17,1)=g;
                        c_array(num).bi(17,2)=0;
                        c_array(num).bi(17,3)=g;
                        
                        c_array(num).bi(18,1)=0;
                        c_array(num).bi(18,2)=g;
                        c_array(num).bi(18,3)=-g;
                        
                        c_array(num).bi(19,1)=g;
                        c_array(num).bi(19,2)=g;
                        c_array(num).bi(19,3)=0;
                        
                        c_array(num).bi(20,1)=g;
                        c_array(num).bi(20,2)=0;
                        c_array(num).bi(20,3)=-g;
                        
                        c_array(num).bi(21,1)=0;
                        c_array(num).bi(21,2)=g;
                        c_array(num).bi(21,3)=g;
                        
                        c_array(num).bi(22,1)=g;
                        c_array(num).bi(22,2)=-g;
                        c_array(num).bi(22,3)=0;
                        
                        c_array(num).bi(23,1)=g;
                        c_array(num).bi(23,2)=0;
                        c_array(num).bi(23,3)=g;
                        
                        c_array(num).bi(24,1)=0;
                        c_array(num).bi(24,2)=g;
                        c_array(num).bi(24,3)=g;

                        
                        f = 1/sqrt(3.0D0);
                        %edge
        %-------------------{1 1 0}<1 -1 0>------------------------
                        c_array(num).ni(1,1)=g;
                        c_array(num).ni(1,2)=g;
                        c_array(num).ni(1,3)=0;
                        
                        c_array(num).ni(2,1)=-g;
                        c_array(num).ni(2,2)=g;
                        c_array(num).ni(2,3)=0;
                        
                        c_array(num).ni(3,1)=g;
                        c_array(num).ni(3,2)=0;
                        c_array(num).ni(3,3)=g;
                        
                        c_array(num).ni(4,1)=g;
                        c_array(num).ni(4,2)=0;
                        c_array(num).ni(4,3)=-g;
                        
                        c_array(num).ni(5,1)=0;
                        c_array(num).ni(5,2)=g;
                        c_array(num).ni(5,3)=g;
                        
                        c_array(num).ni(6,1)=0;
                        c_array(num).ni(6,2)=g;
                        c_array(num).ni(6,3)=-g;
        %-------------------{1 0 0}<0  1 1>------------------------
                        c_array(num).ni(7,1)=1;
                        c_array(num).ni(7,2)=0;
                        c_array(num).ni(7,3)=0;
                        
                        c_array(num).ni(8,1)=1;
                        c_array(num).ni(8,2)=0;
                        c_array(num).ni(8,3)=0;
                        
                        c_array(num).ni(9,1)=0;
                        c_array(num).ni(9,2)=1;
                        c_array(num).ni(9,3)=0;
                        
                        c_array(num).ni(10,1)=0;
                        c_array(num).ni(10,2)=1;
                        c_array(num).ni(10,3)=0;
                        
                        c_array(num).ni(11,1)=0;
                        c_array(num).ni(11,2)=0;
                        c_array(num).ni(11,3)=1;
                        
                        c_array(num).ni(12,1)=0;
                        c_array(num).ni(12,2)=0;
                        c_array(num).ni(12,3)=1;
        %-------------------{1 1 1}<1 -1 0>------------------------
                        c_array(num).ni(13,1)=f;
                        c_array(num).ni(13,2)=f;
                        c_array(num).ni(13,3)=f;
                        
                        c_array(num).ni(14,1)=-f;
                        c_array(num).ni(14,2)=-f;
                        c_array(num).ni(14,3)=-f;
                        
                        c_array(num).ni(15,1)=f;
                        c_array(num).ni(15,2)=f;
                        c_array(num).ni(15,3)=f;
                        
                        c_array(num).ni(16,1)=f;
                        c_array(num).ni(16,2)=-f;
                        c_array(num).ni(16,3)=-f;
                        
                        c_array(num).ni(17,1)=-f;
                        c_array(num).ni(17,2)=f;
                        c_array(num).ni(17,3)=f;
                        
                        c_array(num).ni(18,1)=f;
                        c_array(num).ni(18,2)=-f;
                        c_array(num).ni(18,3)=-f;
                        
                        c_array(num).ni(19,1)=f;
                        c_array(num).ni(19,2)=-f;
                        c_array(num).ni(19,3)=f;
                        
                        c_array(num).ni(20,1)=f;
                        c_array(num).ni(20,2)=-f;
                        c_array(num).ni(20,3)=f;
                        
                        c_array(num).ni(21,1)=-f;
                        c_array(num).ni(21,2)=f;
                        c_array(num).ni(21,3)=-f;
                        
                        c_array(num).ni(22,1)=-f;
                        c_array(num).ni(22,2)=-f;
                        c_array(num).ni(22,3)=f;
                        
                        c_array(num).ni(23,1)=-f;
                        c_array(num).ni(23,2)=-f;
                        c_array(num).ni(23,3)=f;
                        
                        c_array(num).ni(24,1)=f;
                        c_array(num).ni(24,2)=f;
                        c_array(num).ni(24,3)=-f;
%--------------------------------Ran Halite End----------------------------

                  else
                        error('')
                  end %if
% c
% c               
                      e = c_array(num).e;
                      v = c_array(num).nu;
                      u = c_array(num).mu;
                      
                  if (c_array(num).elastic_type == 1) %then
                        c_array(num).elast_flex(1,1)=1/e;
                        c_array(num).elast_flex(1,2)=-v/e;
                        c_array(num).elast_flex(1,3)=-v/e;
                        c_array(num).elast_flex(1,4)=0;
                        c_array(num).elast_flex(1,5)=0;
                        c_array(num).elast_flex(1,6)=0;
                        c_array(num).elast_flex(2,1)=-v/e;
                        c_array(num).elast_flex(2,2)=1/e;
                        c_array(num).elast_flex(2,3)=-v/e;
                        c_array(num).elast_flex(2,4)=0;
                        c_array(num).elast_flex(2,5)=0;
                        c_array(num).elast_flex(2,6)=0;
                        c_array(num).elast_flex(3,1)=-v/e;
                        c_array(num).elast_flex(3,2)=-v/e;
                        c_array(num).elast_flex(3,3)=1/e;
                        c_array(num).elast_flex(3,4)=0;
                        c_array(num).elast_flex(3,5)=0;
                        c_array(num).elast_flex(3,6)=0;
                        c_array(num).elast_flex(4,1)=0;
                        c_array(num).elast_flex(4,2)=0;
                        c_array(num).elast_flex(4,3)=0;
                        c_array(num).elast_flex(4,4)=2*(1+v)/e;
                        c_array(num).elast_flex(4,5)=0;
                        c_array(num).elast_flex(4,6)=0;
                        c_array(num).elast_flex(5,1)=0;
                        c_array(num).elast_flex(5,2)=0;
                        c_array(num).elast_flex(5,3)=0;
                        c_array(num).elast_flex(5,4)=0;
                        c_array(num).elast_flex(5,5)=2*(1+v)/e;
                        c_array(num).elast_flex(5,6)=0;
                        c_array(num).elast_flex(6,1)=0;
                        c_array(num).elast_flex(6,2)=0;
                        c_array(num).elast_flex(6,3)=0;
                        c_array(num).elast_flex(6,4)=0;
                        c_array(num).elast_flex(6,5)=0;
                        c_array(num).elast_flex(6,6)=2*(1+v)/e;
                  elseif (c_array(num).elastic_type == 2) %then
                        c_array(num).elast_flex(1,1)=1/e;
                        c_array(num).elast_flex(1,2)=-v/e;
                        c_array(num).elast_flex(1,3)=-v/e;
                        c_array(num).elast_flex(1,4)=0;
                        c_array(num).elast_flex(1,5)=0;
                        c_array(num).elast_flex(1,6)=0;
                        c_array(num).elast_flex(2,1)=-v/e;
                        c_array(num).elast_flex(2,2)=1/e;
                        c_array(num).elast_flex(2,3)=-v/e;
                        c_array(num).elast_flex(2,4)=0;
                        c_array(num).elast_flex(2,5)=0;
                        c_array(num).elast_flex(2,6)=0;
                        c_array(num).elast_flex(3,1)=-v/e;
                        c_array(num).elast_flex(3,2)=-v/e;
                        c_array(num).elast_flex(3,3)=1/e;
                        c_array(num).elast_flex(3,4)=0;
                        c_array(num).elast_flex(3,5)=0;
                        c_array(num).elast_flex(3,6)=0;
                        c_array(num).elast_flex(4,1)=0;
                        c_array(num).elast_flex(4,2)=0;
                        c_array(num).elast_flex(4,3)=0;
                        c_array(num).elast_flex(4,4)=1/u;
                        c_array(num).elast_flex(4,5)=0;
                        c_array(num).elast_flex(4,6)=0;
                        c_array(num).elast_flex(5,1)=0;
                        c_array(num).elast_flex(5,2)=0;
                        c_array(num).elast_flex(5,3)=0;
                        c_array(num).elast_flex(5,4)=0;
                        c_array(num).elast_flex(5,5)=1/u;
                        c_array(num).elast_flex(5,6)=0;
                        c_array(num).elast_flex(6,1)=0;
                        c_array(num).elast_flex(6,2)=0;
                        c_array(num).elast_flex(6,3)=0;
                        c_array(num).elast_flex(6,4)=0;
                        c_array(num).elast_flex(6,5)=0;
                        c_array(num).elast_flex(6,6)=1/u;
                  elseif (c_array(num).elastic_type == 3) %then
                      % do nothing here
                  else
                      error('')
                  end %if
% c                       Generate (one time) the inverse of the flexibility
                  if (c_array(num).elastic_type ~= 3)
                  c_array(num).elast_stiff = inv(c_array(num).elast_flex);
                  else
                        C11 = c_array(num).mu_o;
                        C33 = c_array(num).D_o;
                        C12 = c_array(num).t_o;
                        C13 = c_array(num).e;
                        C44 = c_array(num).nu;
                        C55 = c_array(num).voche_m;
                      C22 = C11;
                      C23 = C13;   
                      C66 = C55;                    
                        c_array(num).elast_stiff(1,1)=C11;
                        c_array(num).elast_stiff(1,2)=C12;
                        c_array(num).elast_stiff(1,3)=C13;
                        c_array(num).elast_stiff(1,4)=0.d0;
                        c_array(num).elast_stiff(1,5)=0.d0;
                        c_array(num).elast_stiff(1,6)=0.d0;
                        c_array(num).elast_stiff(2,1)=C12;
                        c_array(num).elast_stiff(2,2)=C22;
                        c_array(num).elast_stiff(2,3)=C23;
                        c_array(num).elast_stiff(2,4)=0.d0;
                        c_array(num).elast_stiff(2,5)=0.d0;
                        c_array(num).elast_stiff(2,6)=0.d0;
                        c_array(num).elast_stiff(3,1)=C13;
                        c_array(num).elast_stiff(3,2)=C23;
                        c_array(num).elast_stiff(3,3)=C33;
                        c_array(num).elast_stiff(3,4)=0.d0;
                        c_array(num).elast_stiff(3,5)=0.d0;
                        c_array(num).elast_stiff(3,6)=0.d0;
                        c_array(num).elast_stiff(4,1)=0.d0;
                        c_array(num).elast_stiff(4,2)=0.d0;
                        c_array(num).elast_stiff(4,3)=0.d0;
                        c_array(num).elast_stiff(4,4)=C44;
                        c_array(num).elast_stiff(4,5)=0.d0;
                        c_array(num).elast_stiff(4,6)=0.d0;
                        c_array(num).elast_stiff(5,1)=0.d0;
                        c_array(num).elast_stiff(5,2)=0.d0;
                        c_array(num).elast_stiff(5,3)=0.d0;
                        c_array(num).elast_stiff(5,4)=0.d0;
                        c_array(num).elast_stiff(5,5)=C55;
                        c_array(num).elast_stiff(5,6)=0.d0;
                        c_array(num).elast_stiff(6,1)=0.d0;
                        c_array(num).elast_stiff(6,2)=0.d0;
                        c_array(num).elast_stiff(6,3)=0.d0;
                        c_array(num).elast_stiff(6,4)=0.d0;
                        c_array(num).elast_stiff(6,5)=0.d0;
                        c_array(num).elast_stiff(6,6)=C66;                      
                  end
%                   call mm10_invsym(c_array(num).elast_stiff,6)
% c
                  c_array(num).valid = 1;
          end
      end % Initialize
   end % methods
end
