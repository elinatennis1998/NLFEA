classdef local_work_def < handle % makes it act like pass-by-reference
    % From include_sig_up
% c
% c              declaration of local arrays to be allocated on the
% c              stack for strain-stress updating. this
% c              enables blocks of elements to be processed
% c              in parallel.
% c
% properties (Constant)
%     mxnod=500000;mxel=500000;mxndel=20;max_threads=64;
%  mxlbel=15;mxndof=3;mxgp=14;nstr=6;mxstep=6000;mxlc=2000;
%  mxmat=15;mxlsz=mxnod/3;mxelmp=30;mxstmp=30; mxstc=10;
%  mxoupr=30;ntrc=10;mxconn=50;mxvl=128;mxnmgp=10;
%  mxblsz=128;mxnmbl=2000;nparam=3;mxelpr=42;vclim=28;
%  mxmtpr=200;mxsepr=11;two16=65536;ndim=3;max_tables=20;
%  mxtim=50;nstrs=9;max_procs=1024;mxndpr=8;max_crystals=1000;
%  max_slip_sys=12; max_uhard=12
%    end
    properties
%         ce_0 = zeros(mxvl,mxecor)
%         ce_n = zeros(mxvl,mxecor)
%         ce_mid = zeros(mxvl,mxecor)
%         ce_n1 = zeros(mxvl,mxecor)
%         trnmte = zeros(          mxvl,mxedof,mxndof)
%         det_j = zeros(             mxvl,mxgp)
%         det_j_mid = zeros(         mxvl,mxgp)
%         weights = zeros(mxgp)
%         nxi = zeros(               mxndel,mxgp)
%         neta = zeros(              mxndel,mxgp)
%         nzeta = zeros(             mxndel,mxgp)
%         gama = zeros(          mxvl,3,3,mxgp)
%         gama_mid = zeros(      mxvl,3,3,mxgp)
%         fn = zeros(              mxvl,3,3)
%         fn1 = zeros(             mxvl,3,3)
%         dfn1 = zeros(	            mxvl)
%         vol_block = zeros(       mxvl,8,3)
%         volume_block = zeros(        mxvl)
%         volume_block_0 = zeros(      mxvl)
%         volume_block_n = zeros(      mxvl)
%         volume_block_n1 = zeros(     mxvl)
        jac = zeros(             1,3,3)
%         b = zeros(               mxvl,mxedof,nstr)
%         ue = zeros(                mxvl,mxedof)
%         due = zeros(               mxvl,mxedof)
%         uenh = zeros(              mxvl,mxedof)
%         uen1 = zeros(              mxvl,mxedof)
%         urcs_blk_n = zeros(      mxvl,nstrs,mxgp)
%         urcs_blk_n1 = zeros(     mxvl,nstrs,mxgp)
%         rot_blk_n1 = zeros(      mxvl,9,mxgp)
%         rtse = zeros(            mxvl,nstr,mxgp)
%         elem_hist1 = zeros(100,100,100)
%         elem_hist = zeros(100,100,100)
%         ddtse = zeros(           mxvl,nstr,mxgp)
%         strain_n = zeros(        mxvl,nstr,mxgp)
%         dtemps_node_blk = zeros(   mxvl,mxndel)
%         temps_ref_node_blk = zeros(mxvl,mxndel)
%         temps_node_blk = zeros(    mxvl,mxndel)
%         temps_node_ref_blk = zeros(mxvl,mxndel)
%         cohes_temp_ref = zeros(      mxvl)
%         cohes_dtemp = zeros(         mxvl)
%         cohes_temp_n = zeros(        mxvl)
%         nu_vec = zeros(              mxvl)
%         beta_vec = zeros(            mxvl)
%         h_vec = zeros(               mxvl)
%         e_vec = zeros(               mxvl)
%         sigyld_vec = zeros(          mxvl)
        alpha_vec = zeros(         1,6)
%         e_vec_n = zeros(             mxvl)
%         nu_vec_n = zeros(            mxvl)
%         gp_sig_0_vec = zeros(        mxvl)
%         gp_sig_0_vec_n = zeros(      mxvl)
%         gp_h_u_vec = zeros(          mxvl)
%         gp_h_u_vec_n = zeros(        mxvl)
%         gp_beta_u_vec = zeros(       mxvl)
%         gp_beta_u_vec_n = zeros(     mxvl)
%         gp_delta_u_vec = zeros(      mxvl)
%         gp_delta_u_vec_n = zeros(    mxvl)
        alpha_vec_n = zeros(       1,6)
%         h_vec_n = zeros(             mxvl)
%         n_power_vec = zeros(         mxvl)
%         f0_vec = zeros(              mxvl)
%         eps_ref_vec = zeros(         mxvl)
%         m_power_vec = zeros(         mxvl)
%         q1_vec = zeros(              mxvl)
%         q2_vec = zeros(              mxvl)
%         q3_vec = zeros(              mxvl)
%         nuc_s_n_vec = zeros(         mxvl)
%         nuc_e_n_vec = zeros(         mxvl)
%         nuc_f_n_vec = zeros(         mxvl)

%         ce_0 = zeros(mxvl,mxecor)
%         ce_n = zeros(mxvl,mxecor)
%         ce_mid = zeros(mxvl,mxecor)
%         ce_n1 = zeros(mxvl,mxecor)
%         trnmte = zeros(          mxvl,mxedof,mxndof)
%         det_j = zeros(             mxvl,mxgp)
%         det_j_mid = zeros(         mxvl,mxgp)
%         weights = zeros(mxgp)
%         nxi = zeros(               mxndel,mxgp)
%         neta = zeros(              mxndel,mxgp)
%         nzeta = zeros(             mxndel,mxgp)
%         gama = zeros(          mxvl,3,3,mxgp)
%         gama_mid = zeros(      mxvl,3,3,mxgp)
%         fn = zeros(              mxvl,3,3)
%         fn1 = zeros(             mxvl,3,3)
%         dfn1 = zeros(	            mxvl)
%         vol_block = zeros(       mxvl,8,3)
%         volume_block = zeros(        mxvl)
%         volume_block_0 = zeros(      mxvl)
%         volume_block_n = zeros(      mxvl)
%         volume_block_n1 = zeros(     mxvl)
%         jac = zeros(             mxvl,3,3)
%         b = zeros(               mxvl,mxedof,nstr)
%         ue = zeros(                mxvl,mxedof)
%         due = zeros(               mxvl,mxedof)
%         uenh = zeros(              mxvl,mxedof)
%         uen1 = zeros(              mxvl,mxedof)
        urcs_blk_n = zeros(      1,9,27)
        urcs_blk_n1 = zeros(     1,9,27)
        rot_blk_n1 = zeros(      1,9,27)
%         rtse = zeros(            mxvl,nstr,mxgp)
        elem_hist1 = zeros(1,152,27)
        elem_hist = zeros(1,152,27)
%         ddtse = zeros(           mxvl,nstr,mxgp)
%         strain_n = zeros(        mxvl,nstr,mxgp)
%         dtemps_node_blk = zeros(   mxvl,mxndel)
%         temps_ref_node_blk = zeros(mxvl,mxndel)
%         temps_node_blk = zeros(    mxvl,mxndel)
%         temps_node_ref_blk = zeros(mxvl,mxndel)
%         cohes_temp_ref = zeros(      mxvl)
%         cohes_dtemp = zeros(         mxvl)
%         cohes_temp_n = zeros(        mxvl)
%         nu_vec = zeros(              mxvl)
%         beta_vec = zeros(            mxvl)
%         h_vec = zeros(               mxvl)
%         e_vec = zeros(               mxvl)
%         sigyld_vec = zeros(          mxvl)
%         alpha_vec = zeros(         mxvl,6)
%         e_vec_n = zeros(             mxvl)
%         nu_vec_n = zeros(            mxvl)
%         gp_sig_0_vec = zeros(        mxvl)
%         gp_sig_0_vec_n = zeros(      mxvl)
%         gp_h_u_vec = zeros(          mxvl)
%         gp_h_u_vec_n = zeros(        mxvl)
%         gp_beta_u_vec = zeros(       mxvl)
%         gp_beta_u_vec_n = zeros(     mxvl)
%         gp_delta_u_vec = zeros(      mxvl)
%         gp_delta_u_vec_n = zeros(    mxvl)
%         alpha_vec_n = zeros(       mxvl,6)
%         h_vec_n = zeros(             mxvl)
%         n_power_vec = zeros(         mxvl)
%         f0_vec = zeros(              mxvl)
%         eps_ref_vec = zeros(         mxvl)
%         m_power_vec = zeros(         mxvl)
%         q1_vec = zeros(              mxvl)
%         q2_vec = zeros(              mxvl)
%         q3_vec = zeros(              mxvl)
%         nuc_s_n_vec = zeros(         mxvl)
%         nuc_e_n_vec = zeros(         mxvl)
%         nuc_f_n_vec = zeros(         mxvl)
       %real
       dt=0;total_model_time=0;time_n=0;beta_fact=0;
         block_energy=0;eps_bbar=0;block_plastic_work=0;step_scale_fact=0;
         alpha_dmg=0;ls=0;ll=0;lt=0;
       %integer
       felem=0;elem_type=0;matnum=0;int_order=0;mat_type=0;
                  num_enodes=0;num_enode_dof=0;totdof=0;
                  num_int_points=0;span=0;iter=0;step=0;gpn=0;
                  number_points=0;cohes_type=0;curve_set_number=0;
                  surface=0;hist_size_for_blk=0;iout=0;blk=0;
                  umat_stress_type=0;cep_sym_size=0;num_threads=0;
                  inter_mat=0;macro_sz=0;cp_sz =0;
       %logical
       geo_non_flg=1;bbar_flg=0;trn_e_block=0;
%                   trn_e_flags(mxvl)=0;
                  first=0;material_cut_step=0;signal_flag=0;
%                   adaptive_flag=0;temperatures=0;lnelas_vec(mxvl)=0;
%                   nuc_vec(mxvl)=0;nonlinear_flag(mxvl)=0;
                  allow_cut=0;
                  segmental=0;power_law=0;temps_node_to_process=0;
                  temperatures_ref=0;fgm_enode_props=0;is_cohes_elem=0;
                  linear_displ_elem=0;adjust_const_elem=0;
%                   is_axisymm_elem=0;killed_status_vec(mxvl)=0;
                  block_killed=0;is_umat=0;is_solid_matl=0;is_crys_pls=0;
                  compute_f_bar=0;compute_f_n=0;is_cohes_nonlocal=0;
                  is_inter_dmg=0;
% c     Added stuff for CP
      debug_flag = zeros(1,1);                   % mxvl
      local_tol = zeros(1,1);           % mxvl
      ncrystals = zeros(1,1) ;                   % mxvl
      angle_type = zeros(1,1) ;                  % mxvl
      angle_convention = zeros(1,1);             % mxvl
      c_props = crystal_props;   % mxvl,max_crystals
    end
   
%    methods
%       function LW = local_work_def(step,c_props,w_props) % Initialization method
%           if nargin > 0
%               LW.step  = step;
%               LW.c_props = c_props;
%               LW.angle_type = w_props.angle_type;
%               LW.angle_convention = w_props.angle_convention;
%               LW.angle_type = w_props.angle_type;
%           end
%       end % Initialize
%    end % methods
end % classdef