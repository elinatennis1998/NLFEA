% c
% c     ****************************************************************
% c     *                                                              *
% c     *                 subroutine mm10                              *
% c     *                                                              *
% c     *                       written by : mcm                       *
% c     *                                                              *
% c     *                   last modified: 12/11/13                     *
% c     *                                                              *
% c     *    crystal plasticity stress-strain update                   *
% c     *                                                              *
% c     *                                                              *
% c     ****************************************************************
% c
function [history_np1,history_n, ...
    local_work, gp_temps, gp_temp_inc] = mm10(gp, span,...
    ncrystals, hist_sz, history_n, local_work, uddt,...
    gp_temps, gp_temp_inc, iout,subcycle)
%             use segmental_curves, only: max_seg_points
%             use mm10_defs
[max_slip_sys,max_uhard] = maxparamsCP;
%             implicit integer (a-z)
% $add include_sig_up
%             double precision, intent(in) :: uddt(mxvl,nstr)
%             double precision, intent(in) :: gp_temps(mxvl),
%      &                                      gp_temp_inc(mxvl)
%             integer, intent(in) :: iout
%             integer :: span, ncrystals(mxvl), hist_sz, gp
%             double precision :: history_n(span,hist_sz)
history_np1 = zeros(span,hist_sz);
% c     
%             logical :: debug
%             integer :: i,c,co,cn
%             double precision, dimension(6) :: sig_avg
%             double precision, dimension(6,6) :: tang_avg
%             double precision, dimension(max_slip_sys) :: slip_avg
%             double precision :: t_work_inc, p_work_inc,p_strain_inc
%             type(crystal_props) :: cc_props
%             type(crystal_state) :: cc_n, cc_np1
% c
%            debug = 0 % assuming false = 0
% c
%             if (debug) 
%                 write (*,*) "In mm10"
%             end
% c
% c                 Loop on gauss points
%             local_work.material_cut_step = 0; % assuming .false. = 0
            for i = 1:span
% c                       Initialize the element history if it's step
%                   if (local_work.step == 1) %% NOTE: Tim now handles this
%                   in case isw=40
% %                         if (debug) 
% %                             write(*,*) "Init GP history"
% %                         end
%                         [history_n(i,1:72)] = ...
%                             mm10_init_general_hist(history_n(i,1:72));
%                         [history_n(i,73:75)] = ...
%                             mm10_init_uout_hist(history_n(i,73:75));
%                         [history_n(i,76+12:76+12+max_slip_sys-1)] = ...
%                             mm10_init_slip_hist(history_n(i,76+12:76+12+...
%                             max_slip_sys-1));
%                   end
% c                       Fix a problem with the rotations
% c                       I feel this really should be done elsewhere
                  if (~ local_work.geo_non_flg)
                     local_work.rot_blk_n1(i, 1:9, gp) = ...
                         ([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
                  end
% 
%                   if (debug)
%                       write (*,*) "Updating local ", i
%                   end
                  sig_avg = 0.0;
                  slip_avg = 0.0;
                  t_work_inc = 0.0;
                  p_work_inc = 0.0;
                  p_strain_inc = 0.0;
                  tang_avg = 0.0;
% c                       Loop on crystals
                  for c = 1:ncrystals(i)
                      co = 76+12+max_slip_sys+(c-1)*(30+max_slip_sys+3*max_uhard);
                      cn = 76+12+max_slip_sys+(c)*(30+max_slip_sys+3*max_uhard)-1;
%                       if (debug)
%                           write(*,*) "Setting up properties"
%                       end
                      cc_props...
                          = mm10_init_cc_props(local_work.c_props(i,c)...
                          ,local_work.angle_type(i),...
                          local_work.angle_convention(i));
                      cc_props.out = iout;
% c
%                       if (local_work.step == 1) %% NOTE: Tim now handles this
% %                   in case isw=40
% %                         if (debug) 
% %                             write(*,*) "Init history 0"
% %                         end
%                         [cc_props,local_work.c_props(i,c).init_angles(1:3)...
%                             ,history_n(i,co:cn)] = ...
%                             mm10_init_cc_hist0(cc_props,...
%                             local_work.c_props(i,c).init_angles(1:3)...
%                             ,history_n(i,co:cn));
%                       end
%                       if (debug) 
%                           write(*,*) "Copying n to struct"
%                       end
                      [history_n(i,co:cn),history_n(i,37:63),...
                          history_n(i,64:72),cc_props, cc_n] = ...
                          mm10_copy_cc_hist(history_n(i,co:cn),...
                          history_n(i,37:63),history_n(i,64:72)...
                          ,cc_props,history_n(i,76:76+11));
% c
%                       [local_work.rot_blk_n1(i,1:9,gp), uddt(i,1:6), ...
%                           local_work.dt, gp_temps(i), ...
%                           local_work.step,i+local_work.felem, gp, ...
                          cc_np1 = mm10_setup_np1...
                          (local_work.rot_blk_n1(i,1:9,gp), ...
                          uddt(i,1:6), local_work.dt, gp_temps(i), ...
                          local_work.step,i+local_work.felem, gp);
% c
%                       if (debug) 
%                           write (*,*) "Updating crystal ", c
%                       end
                      [cc_props, cc_np1, cc_n, ...
                          local_work.material_cut_step] = ...
                          mm10_solve_crystal(cc_props, cc_np1, cc_n, ...
                          local_work.material_cut_step, iout, 0,subcycle); % .false. replaced by 0 in last variable in argument
% c
                      if (local_work.material_cut_step)
                          return
                      end
% c                       
% c                     Add stuff into the average
                      sig_avg = sig_avg + cc_np1.stress;
                      tang_avg = tang_avg + cc_np1.tangent;
                      slip_avg = slip_avg + cc_np1.slip_incs;
                      t_work_inc = t_work_inc + cc_np1.work_inc;
                      p_work_inc = p_work_inc + cc_np1.p_work_inc;
                      p_strain_inc = p_strain_inc + cc_np1.p_strain_inc;
% c
% c                     Store the CP history for this crystal
                      history_np1(i,co:cn) ...
                          = mm10_store_cryhist(cc_props, cc_np1, ...
                          cc_n, history_np1(i,co:cn));
                  end
% c
% c                 Do the division for the averages
                  sig_avg = sig_avg / (ncrystals(i));
                  tang_avg = tang_avg / (ncrystals(i));
                  slip_avg = slip_avg / (ncrystals(i));
                  t_work_inc = t_work_inc / (ncrystals(i));
                  p_work_inc = p_work_inc / (ncrystals(i));
                  p_strain_inc = p_strain_inc / (ncrystals(i));
% c
% c                 Actually store the GP data
% c
                  [local_work.urcs_blk_n1(i,1:6,gp), ...
                      history_np1(i,1:36), ...
                      history_np1(i,76+12:76+12+max_slip_sys-1), ...
                      local_work.urcs_blk_n1(i,7:9,gp), ...
                      history_np1(i,64:72)] = mm10_store_gp(sig_avg, ...
                      tang_avg, ...
                      slip_avg, ...
                      history_n(i,76+12:76+12+max_slip_sys-1), ...
                      t_work_inc, ...
                      p_work_inc, p_strain_inc, ...
                      local_work.urcs_blk_n(i,7:9,gp), ...
                      local_work.rot_blk_n1(i,1:9,gp));
            end
% c
%             return
end           
                  
                  
                  
                  
                  
                  
                  
                  