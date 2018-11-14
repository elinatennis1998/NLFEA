% c
% c     ****************************************************************
% c     *                                                              *
% c     *                   subroutine setup_mm10_rknstr               *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified : 03/23/2012                 *
% c     *                                                              *
% c     *     set up material model #10 (crystal plasticity)           *
% c     *     for stress updating: values constant across all g. pts.  *
% c     *                                                              *
% c     ****************************************************************
% c
% c
function local_work = setup_mm10_rknstr( span,c_array, props, elem, dt)
%       use segmental_curves
%       use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
%       use crystal_data, only : c_array, angle_input, crystal_input,
%      &                              data_offset
% c
%       implicit integer (a-z)
% $add param_def
% c
% c                    parameter declarations
% c
%       real    props(mxelpr,mxvl)
%       logical lprops(mxelpr,mxvl)
%       integer iprops(mxelpr,mxvl)
% $add include_sig_up
% c
% c                    local
% c
%       logical adaptive
%       integer :: ctotal, c, cnum, s
%       double precision, dimension(3) :: angles
%       character :: aconv*5
%       character :: atype*7
%       double precision, dimension(3) :: bs, ns
%       double precision, dimension(3,3) :: A
%       double precision, dimension(6,6) :: Rstiff, temp
%       integer :: elnum, osn
% c
      ctotal = 0;
      local_work = local_work_def;
% c
%       matnum = local_work.matnum;
% c
      for i = 1: span
% c
%            local_work.alpha_vec(i,1) = props(9,i);
%            local_work.alpha_vec(i,2) = props(13,i);
%            local_work.alpha_vec(i,3) = props(34,i);
%            local_work.alpha_vec(i,4) = props(35,i);
%            local_work.alpha_vec(i,5) = props(36,i);
%            local_work.alpha_vec(i,6) = props(37,i);
% % c
%            local_work.alpha_vec_n(i,1) = props(9,i);
%            local_work.alpha_vec_n(i,2) = props(13,i);
%            local_work.alpha_vec_n(i,3) = props(34,i);
%            local_work.alpha_vec_n(i,4) = props(35,i);
%            local_work.alpha_vec_n(i,5) = props(36,i);
%            local_work.alpha_vec_n(i,6) = props(37,i);
% c
%            local_work.debug_flag(i) = lmtprp(13,matnum);
%            local_work.local_tol(i) = dmatprp(100,matnum);
% c                 May eventually change this to allow for different # of
% c                 crystals in block
           local_work.ncrystals(i) = props.ncrystals;
% c
           local_work.angle_convention(i) = props.angle_convention;
           local_work.angle_type(i) = props.angle_type;
% c

% c                 VERY IMPORTANT LOOP
% c                 Will need to (in the near future) extract crystal
% c                 and orientation information.  Also possibly change
% c                 to allow for variable number of crystals in each element
% c                 block.
% c
%            elnum = local_work.felem+i-1;
           local_work.felem = elem - 1;
           local_work.dt = dt;
% c
% c
           for c=1:local_work.ncrystals(i)
% c                       Get the local crystal number
%                   if (imatprp(104,matnum)  ==  1) %then
                        cnum = 1;%imatprp(105,matnum);
%                   elseif (imatprp(104,matnum)  ==  2) %then
%                         osn = data_offset(elnum);
%                         cnum = crystal_input(osn,c);
% % c                       Couldn't do this earlier, so check here
%                         if ((cnum .gt. max_crystals) || ...
%                               (cnum < 0)) %then
%                          error('("Crystal ", i3, " not valid")')
% %                               cnum
% %                               call die_gracefully
%                         end %if
%                   else
%                         error('')
%                         call die_gracefully
%                   end %if

% c                       Get the local orientation
%                   if (imatprp(107,matnum)  ==  1) %then
                        angles(1) =  props.angles(1);
                        angles(2) =  props.angles(2);
                        angles(3) =  props.angles(3);
%                   elseif (imatprp(107,matnum)  ==  2) %then
%                         osn = data_offset(elnum);
%                         angles(1:3) = angle_input(osn,c,1:3);
%                   else
%                         error('')
%                         call die_gracefully
%                   end %if
% c                       Now we have the properties, we just need to extract
% c                       into our local structure
                  local_work.c_props(i,c).stiffness = ...
                        c_array(cnum).elast_stiff;
                  local_work.c_props(i,c).init_angles = angles;
                  local_work.c_props(i,c).nslip = c_array(cnum).nslip;
                  local_work.c_props(i,c).rate_n = c_array(cnum).harden_n;
                  local_work.c_props(i,c).tau_hat_y ...
                        = c_array(cnum).tau_hat_y;
                  local_work.c_props(i,c).G_0_y = c_array(cnum).g_o_y;
                  local_work.c_props(i,c).tau_hat_v ...
                        = c_array(cnum).tau_hat_v;
                  local_work.c_props(i,c).G_0_v = c_array(cnum).g_o_v;
                  local_work.c_props(i,c).tau_a = c_array(cnum).tau_a;
                  local_work.c_props(i,c).burgers = c_array(cnum).b;
                  local_work.c_props(i,c).p_v = c_array(cnum).p_v;
                  local_work.c_props(i,c).q_v = c_array(cnum).q_v;
                  local_work.c_props(i,c).p_y = c_array(cnum).p_y;
                  local_work.c_props(i,c).q_y = c_array(cnum).q_y;
                  local_work.c_props(i,c).boltzman = c_array(cnum).boltz;
                  local_work.c_props(i,c).theta_0 = ...
                        c_array(cnum).theta_o;
                  local_work.c_props(i,c).eps_dot_0_v = ...
                        c_array(cnum).eps_dot_o_v;
                  local_work.c_props(i,c).eps_dot_0_y = ...
                        c_array(cnum).eps_dot_o_y;
                  local_work.c_props(i,c).mu_0 = ...
                        c_array(cnum).mu_o;
                  local_work.c_props(i,c).D_0  = ...
                        c_array(cnum).D_o;
                  local_work.c_props(i,c).T_0 = c_array(cnum).t_o;
                  local_work.c_props(i,c).tau_a = c_array(cnum).tau_a;
                  local_work.c_props(i,c).k_0 = c_array(cnum).k_o;
                  local_work.c_props(i,c).h_type = ...
                        c_array(cnum).h_type;
                  local_work.c_props(i,c).s_type = ...
                        c_array(cnum).slip_type;
                  local_work.c_props(i,c).u1 = c_array(cnum).u1;
                  local_work.c_props(i,c).u2 = c_array(cnum).u2;
                  local_work.c_props(i,c).u3 = c_array(cnum).u3;
                  local_work.c_props(i,c).u4 = c_array(cnum).u4;
                  local_work.c_props(i,c).u5 = c_array(cnum).u5;
                  local_work.c_props(i,c).u6 = c_array(cnum).u6;
                  local_work.c_props(i,c).tau_y = c_array(cnum).tau_y;
                  local_work.c_props(i,c).tau_v = c_array(cnum).tau_v;
                  local_work.c_props(i,c).voche_m =  ...
                        c_array(cnum).voche_m;
                    % New CP model parameters
                  local_work.c_props(i,c).c1 = c_array(cnum).c1;
                  local_work.c_props(i,c).c2 = c_array(cnum).c2;
                  local_work.c_props(i,c).c3 = c_array(cnum).c3;
                  local_work.c_props(i,c).c4 = c_array(cnum).c4;
                  local_work.c_props(i,c).c5 = c_array(cnum).c5;
                  local_work.c_props(i,c).c6 = c_array(cnum).c6;
                  local_work.c_props(i,c).c7 = c_array(cnum).c7;
                  local_work.c_props(i,c).c8 = c_array(cnum).c8;
                  local_work.c_props(i,c).mfp = c_array(cnum).mfp;
                  local_work.c_props(i,c).v0 = c_array(cnum).v0;
                  local_work.c_props(i,c).Qslip = c_array(cnum).Qslip;
                  local_work.c_props(i,c).Qbulk = c_array(cnum).Qbulk;
                  local_work.c_props(i,c).bi = c_array(cnum).bi;
                  local_work.c_props(i,c).ni = c_array(cnum).ni;
                  local_work.c_props(i,c).usefsolve ...
                        = c_array(cnum).usefsolve;
                  local_work.c_props(i,c).plugin ...
                        = c_array(cnum).plugin;
                  local_work.c_props(i,c).miter ...
                      = c_array(cnum).miter;
                  
                  % copy generalized material parameters to local work
                  local_work.c_props(i,c).cp_001 = c_array(cnum).cp_001;
                  local_work.c_props(i,c).cp_002 = c_array(cnum).cp_002;
                  local_work.c_props(i,c).cp_003 = c_array(cnum).cp_003;
                  local_work.c_props(i,c).cp_004 = c_array(cnum).cp_004;
                  local_work.c_props(i,c).cp_005 = c_array(cnum).cp_005;
                  local_work.c_props(i,c).cp_006 = c_array(cnum).cp_006;
                  local_work.c_props(i,c).cp_007 = c_array(cnum).cp_007;
                  local_work.c_props(i,c).cp_008 = c_array(cnum).cp_008;
                  local_work.c_props(i,c).cp_009 = c_array(cnum).cp_009;
                  local_work.c_props(i,c).cp_010 = c_array(cnum).cp_010;
                  local_work.c_props(i,c).cp_011 = c_array(cnum).cp_011;
                  local_work.c_props(i,c).cp_012 = c_array(cnum).cp_012;
                  local_work.c_props(i,c).cp_013 = c_array(cnum).cp_013;
                  local_work.c_props(i,c).cp_014 = c_array(cnum).cp_014;
                  local_work.c_props(i,c).cp_015 = c_array(cnum).cp_015;
                  local_work.c_props(i,c).cp_016 = c_array(cnum).cp_016;
                  local_work.c_props(i,c).cp_017 = c_array(cnum).cp_017;
                  local_work.c_props(i,c).cp_018 = c_array(cnum).cp_018;
                  local_work.c_props(i,c).cp_019 = c_array(cnum).cp_019;
                  local_work.c_props(i,c).cp_020 = c_array(cnum).cp_020;
                  local_work.c_props(i,c).cp_021 = c_array(cnum).cp_021;
                  local_work.c_props(i,c).cp_022 = c_array(cnum).cp_022;
                  local_work.c_props(i,c).cp_023 = c_array(cnum).cp_023;
                  local_work.c_props(i,c).cp_024 = c_array(cnum).cp_024;
                  local_work.c_props(i,c).cp_025 = c_array(cnum).cp_025;
                  local_work.c_props(i,c).cp_026 = c_array(cnum).cp_026;
                  local_work.c_props(i,c).cp_027 = c_array(cnum).cp_027;
                  local_work.c_props(i,c).cp_028 = c_array(cnum).cp_028;
                  local_work.c_props(i,c).cp_029 = c_array(cnum).cp_029;
                  local_work.c_props(i,c).cp_030 = c_array(cnum).cp_030;
                  local_work.c_props(i,c).cp_031 = c_array(cnum).cp_031;
                  local_work.c_props(i,c).cp_032 = c_array(cnum).cp_032;
                  local_work.c_props(i,c).cp_033 = c_array(cnum).cp_033;
                  local_work.c_props(i,c).cp_034 = c_array(cnum).cp_034;
                  local_work.c_props(i,c).cp_035 = c_array(cnum).cp_035;
                  local_work.c_props(i,c).cp_036 = c_array(cnum).cp_036;
                  local_work.c_props(i,c).cp_037 = c_array(cnum).cp_037;
                  local_work.c_props(i,c).cp_038 = c_array(cnum).cp_038;
                  local_work.c_props(i,c).cp_039 = c_array(cnum).cp_039;
                  local_work.c_props(i,c).cp_040 = c_array(cnum).cp_040;
                  local_work.c_props(i,c).cp_041 = c_array(cnum).cp_041;
                  local_work.c_props(i,c).cp_042 = c_array(cnum).cp_042;
                  local_work.c_props(i,c).cp_043 = c_array(cnum).cp_043;
                  local_work.c_props(i,c).cp_044 = c_array(cnum).cp_044;
                  local_work.c_props(i,c).cp_045 = c_array(cnum).cp_045;
                  local_work.c_props(i,c).cp_046 = c_array(cnum).cp_046;
                  local_work.c_props(i,c).cp_047 = c_array(cnum).cp_047;
                  local_work.c_props(i,c).cp_048 = c_array(cnum).cp_048;
                  local_work.c_props(i,c).cp_049 = c_array(cnum).cp_049;
                  local_work.c_props(i,c).cp_050 = c_array(cnum).cp_050;
                  local_work.c_props(i,c).cp_051 = c_array(cnum).cp_051;
                  local_work.c_props(i,c).cp_052 = c_array(cnum).cp_052;
                  local_work.c_props(i,c).cp_053 = c_array(cnum).cp_053;
                  local_work.c_props(i,c).cp_054 = c_array(cnum).cp_054;
                  local_work.c_props(i,c).cp_055 = c_array(cnum).cp_055;
                  local_work.c_props(i,c).cp_056 = c_array(cnum).cp_056;
                  local_work.c_props(i,c).cp_057 = c_array(cnum).cp_057;
                  local_work.c_props(i,c).cp_058 = c_array(cnum).cp_058;
                  local_work.c_props(i,c).cp_059 = c_array(cnum).cp_059;
                  local_work.c_props(i,c).cp_060 = c_array(cnum).cp_060;
                  local_work.c_props(i,c).cp_061 = c_array(cnum).cp_061;
                  local_work.c_props(i,c).cp_062 = c_array(cnum).cp_062;
                  local_work.c_props(i,c).cp_063 = c_array(cnum).cp_063;
                  local_work.c_props(i,c).cp_064 = c_array(cnum).cp_064;
                  local_work.c_props(i,c).cp_065 = c_array(cnum).cp_065;
                  local_work.c_props(i,c).cp_066 = c_array(cnum).cp_066;
                  local_work.c_props(i,c).cp_067 = c_array(cnum).cp_067;
                  local_work.c_props(i,c).cp_068 = c_array(cnum).cp_068;
                  local_work.c_props(i,c).cp_069 = c_array(cnum).cp_069;
                  local_work.c_props(i,c).cp_070 = c_array(cnum).cp_070;
                  local_work.c_props(i,c).cp_071 = c_array(cnum).cp_071;
                  local_work.c_props(i,c).cp_072 = c_array(cnum).cp_072;
                  local_work.c_props(i,c).cp_073 = c_array(cnum).cp_073;
                  local_work.c_props(i,c).cp_074 = c_array(cnum).cp_074;
                  local_work.c_props(i,c).cp_075 = c_array(cnum).cp_075;
                  local_work.c_props(i,c).cp_076 = c_array(cnum).cp_076;
                  local_work.c_props(i,c).cp_077 = c_array(cnum).cp_077;
                  local_work.c_props(i,c).cp_078 = c_array(cnum).cp_078;
                  local_work.c_props(i,c).cp_079 = c_array(cnum).cp_079;
                  local_work.c_props(i,c).cp_080 = c_array(cnum).cp_080;
                  local_work.c_props(i,c).cp_081 = c_array(cnum).cp_081;
                  local_work.c_props(i,c).cp_082 = c_array(cnum).cp_082;
                  local_work.c_props(i,c).cp_083 = c_array(cnum).cp_083;
                  local_work.c_props(i,c).cp_084 = c_array(cnum).cp_084;
                  local_work.c_props(i,c).cp_085 = c_array(cnum).cp_085;
                  local_work.c_props(i,c).cp_086 = c_array(cnum).cp_086;
                  local_work.c_props(i,c).cp_087 = c_array(cnum).cp_087;
                  local_work.c_props(i,c).cp_088 = c_array(cnum).cp_088;
                  local_work.c_props(i,c).cp_089 = c_array(cnum).cp_089;
                  local_work.c_props(i,c).cp_090 = c_array(cnum).cp_090;
                  local_work.c_props(i,c).cp_091 = c_array(cnum).cp_091;
                  local_work.c_props(i,c).cp_092 = c_array(cnum).cp_092;
                  local_work.c_props(i,c).cp_093 = c_array(cnum).cp_093;
                  local_work.c_props(i,c).cp_094 = c_array(cnum).cp_094;
                  local_work.c_props(i,c).cp_095 = c_array(cnum).cp_095;
                  local_work.c_props(i,c).cp_096 = c_array(cnum).cp_096;
                  local_work.c_props(i,c).cp_097 = c_array(cnum).cp_097;
                  local_work.c_props(i,c).cp_098 = c_array(cnum).cp_098;
                  local_work.c_props(i,c).cp_099 = c_array(cnum).cp_099;
                  local_work.c_props(i,c).cp_100 = c_array(cnum).cp_100;
                    
% c                 Call a helper to get the crystal -> reference rotation
                  if (local_work.angle_type(i)  ==  1) %then
                        atype = 'degrees';
                  elseif (local_work.angle_type(i)  ==  2) %then
                        atype = 'radians';
                  else
                        error('')
%                         call die_gracefully
                  end %if

                  if (local_work.angle_convention(i)  ==  1) %then
                        aconv='kocks';
                  else
                        error('')
%                         call die_gracefully
                  end %if
                  
                  
                  switch local_work.c_props(i,c).h_type % Allocate number of hardening variables
                      case 1 %Voche
                          local_work.c_props(i,c).num_hard = 1;
                      case 2 %MTS
                          local_work.c_props(i,c).num_hard = 1;
                      case 3 %User
                          local_work.c_props(i,c).num_hard = 0;
                      case 4 %ORNL
                          local_work.c_props(i,c).num_hard = local_work.c_props(i,c).nslip;
                      case 5 %AADP
                          local_work.c_props(i,c).num_hard = 24;
            case 6 %MTS_omar copy
                local_work.c_props(i,c).num_hard = 1;
                      case 7 %MRR
                          local_work.c_props(i,c).num_hard = local_work.c_props(i,c).nslip;
                      case 8 %Dyson
                          local_work.c_props(i,c).num_hard = 1;
                      case 9 % Ran_HCP
                          local_work.c_props(i,c).num_hard = local_work.c_props(i,c).nslip;
                      case 12 % Plugin library
                          local_work.c_props(i,c).num_hard = mm10_plugin_lib(1,local_work.c_props(i,c),0);
                      case 14 % halite
                          local_work.c_props(i,c).num_hard = local_work.c_props(i,c).nslip;
                  end

                  local_work.c_props(i,c).g = ...
    mm10_rotation_matrix(local_work.c_props(i,c).init_angles, ...
                        aconv, atype);
                    % SOMEHOW, ms and qs get loaded with garbage now inside
                    % of the code, even when they should be initialize to
                    % zero. So, reset them here.
                    local_work.c_props(i,c).ns = zeros(3,1);
                    local_work.c_props(i,c).qs = zeros(3,1);
                    local_work.c_props(i,c).ms = zeros(6,1);
% c
% c                 Now that we have that, we can set up and rotate our
% c                 orientation tensors
                  for s=1:local_work.c_props(i,c).nslip
                        bs = (transpose( ...
                         local_work.c_props(i,c).g)* ...
                         c_array(cnum).bi(s,:)');
                        ns = (transpose( ...
                         local_work.c_props(i,c).g)* ...
                         c_array(cnum).ni(s,:)');
                        local_work.c_props(i,c).ns(1:3,s) = ...
                              ns;
%                         A = spread(bs,dim=2,ncopies=size(ns))*spread( ...
%                               ns,dim=1,ncopies=size(bs))
                        A = (bs*[1 1 1]).*([1 1 1]'*ns');
                        local_work.c_props(i,c).ms(:,s) = mm10_ET2EV(0.5*(A+transpose(A)));
                        local_work.c_props(i,c).qs(:,s) = mm10_WT2WV(0.5*(A-transpose(A)));
                  end %do

% c           We can also rotate forward our stiffness tensor
            Rstiff = mm10_RT2RVE( transpose(local_work.c_props(i,c).g) );
            local_work.c_props(i,c).stiffness = ( ...
                  Rstiff*( ...
                  local_work.c_props(i,c).stiffness* ...
                  transpose(Rstiff)));

           end %do

           ctotal = ctotal + local_work.ncrystals(i);
% c
      end %do

% c
% c                   determine if material model can call for a
% c                   reduction in the adaptive step size
% c
      local_work.allow_cut = 1;
% c
% c
%       return

%  9501 format(/,1x,'>>>> Not implemented yet in rknstr setup!'/)
%  9502 format(/,1x,'>>>> System error: unexpected input type in rknstr!',
%      &            ' Aborting.'/)
%  9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
%      &            ' Aborting.'/)
%  9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
%      &            ' Aborting.'/)
      end