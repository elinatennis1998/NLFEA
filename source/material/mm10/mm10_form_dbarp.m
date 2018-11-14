% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c
% c           Form d_bar_p
function [props, np1, n, vec1, vec2, stress, tt, dbar] = ...
    mm10_form_dbarp(props, np1, n, vec1, vec2, stress, tt)
%      use mm10_defs
%       implicit none
% c
%       type(crystal_props) :: props
%       type(crystal_state) :: np1, n
%       double precision, dimension(6) :: stress
%       double precision, dimension(6) :: dbar
%       double precision :: tt
% c
%       integer :: i
% c
%       double precision :: mm10_slipinc
% c
      if (props.h_type == 1) % voche
        dbarmat = np1.dg * (ones(6,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.ms;
        dbar = sum(dbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 2) % MTS
      rss = stress*np1.ms;
      rss2 = rss.*sign(real(rss));
        dbarmat = np1.dg * (ones(6,1)*((rss2/tt).^(props.rate_n).*sign(real(rss)))).*np1.ms;
        dbar = sum(dbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 3) % User
        dbarmat = np1.dg * (ones(6,1)*((stress*np1.ms/tt).^(props.rate_n).*sign(stress*np1.ms))).*np1.ms;
        dbar = sum(dbarmat(:,1:props.nslip),2);
      elseif (props.h_type == 4) % ORNL
        dbar = np1.ms*vec1(1:props.nslip);
      elseif (props.h_type == 5) % AADP
        dbar = zeros(6,1);
        for i = 1:props.nslip
            dbar = dbar + mm10_slipinc_aadp(props, np1, n, stress, tt, i)*...
                  np1.ms(1:6,i);
        end
	  elseif (props.h_type == 6) % MTS_Omar
		rss = stress*np1.ms - np1.backstress_omar;		
		rss2 = rss.*sign(real(rss));
        dbarmat = np1.dg * (ones(6,1)*((rss2/tt).^(props.rate_n).*sign(real(rss)))).*np1.ms;
        dbar = sum(dbarmat(:,1:props.nslip),2);
	  elseif (props.h_type == 7) % MRR
			dbar = np1.ms*vec1(1:props.nslip);
%         dbar = zeros(6,1);
%         gam = vec1*0;
%         for i = 1:props.nslip
%             gam(i) = mm10_slipinc_mrr(props, np1, n, stress, tt, i);
%             dbar = dbar + gam(i)*...
%                   np1.ms(1:6,i);
%         end
%         if norm(dbar2-dbar)>1e-24
%             dbar;
%         else
%             dbar = dbar2;
%         end
      elseif (props.h_type == 8) % DYS
        dbar = zeros(6,1);
        for i = 1:props.nslip
            dbar = dbar + mm10_slipinc_dys(props, np1, n, stress, tt, i)*...
                  np1.ms(1:6,i);
        end
      elseif (props.h_type == 9) % Ran_HCP
        dbar = np1.ms*vec1(1:props.nslip);
      elseif (props.h_type == 12) % Plugin Library
          dbar = mm10_plugin_lib(6,props,{props, np1, n, vec1, vec2, stress, tt});
      elseif (props.h_type == 14) % halite
        dbar = np1.ms*vec1(1:props.nslip);
      else
        [props] = mm10_unknown_hard_error(props);
      end
% c
%       return
end