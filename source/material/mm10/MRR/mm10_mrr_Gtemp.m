% c
% c *****************************************************************************
% c *                                                                           *
% c *         Small constitutive functions                                      *
% c *                                                                           *
% c *****************************************************************************
% c
% c       Compute temperature dependent G modulus
function G = mm10_mrr_Gtemp(theta)

        
        % Compute the shear modulus using Roter's function
        K11=123.323+6.7008e-8*theta^3-1.1342e-4*theta^2-7.8788e-3*theta;
        K12=70.6512+4.4105e-8*theta^3-7.5498e-5*theta^2-3.9992e-3*theta;
        K44=31.2071+7.0477e-9*theta^3-1.2136e-5*theta^2-8.3274e-3*theta;
        G = 1/3*(K11-K12+K44)*1e9;
      
end
