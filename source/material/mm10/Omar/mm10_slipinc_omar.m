% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       from mm10_form.f
% c
% c           Calculate the slip increment along system i  
function slipinc = ...
    mm10_slipinc_omar(rate_n, dg, ms, stress, tt)

      rs = stress*ms;
      slipinc = dg * (rs/tt)^(rate_n)*sign(rs);