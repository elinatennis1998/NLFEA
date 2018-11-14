% Tim Truster
% 11/02/2012
% Time-stepping loop with Newton-Raphson
% Modified 11/02/2012 to use solution algorithm from Wriggers NLFEM for
% arc-length initialization, and using Mark Messner's idea of using the
% previous values of the prescribed displacements combined with
% F~=F-Kdf*(gBC-gBC_n) for the first iteration.

% 03/02/2015 - added adaptive time step cutting using Warp3d approach
% Now: set the MASTER flag failFEelem == 1 inside an element routine if
% user desires to terminate a NR loop and cut the time step size.
% Also has a divergence check for Rnorm.

% modified 1/17/13 to get arc-length to work again for CohShearAxialTest2


if ~exist('tstep','var')
    tstep = 1;
end

%% continue 10
% Loop for time-history of problem solution

while step < stepmax && ~(transient==7 && lamdamax < lamda + 10*eps) %lamda < lamdamax %lamda <= lamdamax 

 % Advance time marching parameters
    
%  lamda = lamda + s_del_lamda_n1;
 step = step + 1

 % Store data from step n-1 in safe place
 if exist('tstep','var')
     tstep_step = tstep;
 end
 ModelDxn_1_step = ModelDxn_1;
%  ModelVxn_1_step = ModelVxn_1;
%  ModelAxn_1_step = ModelAxn_1;
 if exist('ex_del_ModelDx','var')
     ex_del_ModelDx_step = ex_del_ModelDx;
 end
%     gBC_n_step = gBC_n;
%     hrvec_step = hrvec;
    
 % Set up adaptive time stepping
 stepnotdone1 = 1;
 numincrem1 = 1;
 incremstart1 = 1;
 numincrem2 = 1;
    
    
 %% continue 20
 % Loop to obtain convergence of load step, using adaptive step cutting
 while stepnotdone1
    
      
  %% continue 30
  for increment1 = incremstart1:numincrem1 % incrementation for level 1
    
   if printRnorm && numincrem1 > 1
       increment1
   end
       
    
   %% continue 40
   for increment2 = 1:numincrem2 % incrementation for level 2
        
    if printRnorm && numincrem2 > 1
        increment2
    end
       
    % Set up increment flags
    incrementfailed = 0;
    iter = 0;
    initia = 1; %flag to prescribe an elastic predictor step
    % Master element fail flag - set this to TRUE inside an element routine
    failFEelem = 0;
    diverge_flag = 0;
    
    
    
    
    %% switch 50
    % Compute any required initial quantities
    switch transient
        
        %%
        case 0
            
            s_del_ModelDx = zeros(neq,1);
            if numincrem1 > 1 % Interpolate between steps for current adapted increment
              if step > 1
                lamda = mults(step-1)+(mults(step)-mults(step-1))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
              else
                lamda = mults(step)*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
              end
            else
              lamda = mults(step);
            end
            residtol = 1;
%              rarc = 0;
            
            GetAllLoads
            isw = 3;
            FormFE
            Kn = Kdd11;
            GetAllBCs
            Fd1 = Fd1 - Kdf1*(gBC-gBC_n);
            
            del_ModelDx = Kn\Fd1;
            
            initia = 0;
           
            % Perform line search
            Rnorm = norm(Fd1);
            Residtol = 0;
            linesearchFE

            ModelDx = ModelDx + sLS*del_ModelDx;

            initia = 0;
                
            % Form tangent matrix and residual vector

            GetAllLoads
            GetAllBCs;
%               Kn = Kdd11;
            if mnr == 1
                isw = 6;
                FormFE
                if failFEelem % an element failed
                    break % out of initial step
                end
            else
                isw = 3;
                FormFE%LNL
                if failFEelem % an element failed
                    break % out of initial step
                end
                Kn = Kdd11;
%                     detKn = det(Kdd11);
            end
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e | %+1.7e | %+1.7e\n',Rnorm,log(Rnorm)/log(10),0)
                Rnorm1 = Rnorm;
            end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            %compute initial tangent and residual

        %%
        case {1,2}
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            if enercons == 1
                lener = 0;
            end            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            % User loads
            GetAllLoads
            isw = 3;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                if numbernonlinear == 0 && exist('factorize','file')
%                 yvec = (KnR')\(Schol'*Fd1);
%                 xvec = KnR\yvec;
%                 del_ModelDx = Schol*xvec;
                del_ModelDx = KnR\Fd1;
                elseif numbernonlinear == 0
                del_ModelDx = Kdd11\Fd1;
                else
%                 Kn = Kprime + coeff5*Kdd11;
                del_ModelDx = Kdd11\Fd1;
                end

                % Update field values with increments

                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;

                % Perform line search on increments

                % Form tangent matrix and residual vector

                GetAllLoads
    %           Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            
        %%
        case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
            
            lamda = lamda + s_del_lamda_n1;
            
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            if enercons == 1 % initialize LM for energy to zero
                lener = 0;
            end
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            % first initialize assuming lener=0, meaning linear elasticity

            GetAllLoads
            isw = 3;
            FormFE
            Kn = Kprime + coeffk*Kdd11;
                    
            N = -(Fd1 -(coeffkl2*Fext1 - F1n_1) + Mdd11*ModelAx); %flip sign, subtract off other stuff done in element routine
            A11 = (4/tstep^2*Mdd11 + Kdd11);%(1+lener)*(4/tstep^2*Mdd11 + Kdd11);
            b1 = -N;%-(1+lener)*(4/tstep^2*Mdd11*ModelDx + N) + Mdd11*(ModelAxn_1 + 4/tstep^2*((1+lener)*ModelDxn_1 + tstep*(1+lener/2)*ModelVxn_1));
            y = A11\b1;
            del_ModelDx = y;
            ModelDx = ModelDx + coeffd*del_ModelDx;
            ModelVx = ModelVx + coeffv*del_ModelDx;
            ModelAx = ModelAx + coeffa*del_ModelDx;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Kn = Kprime + coeffk*Kdd11;
                    
            N = -(Fd1 -(coeffkl2*Fext1 - F1n_1) + Mdd11*ModelAx); %flip sign, subtract off other stuff done in element routine
            A11 = (1+lener)*(4/tstep^2*Mdd11 + Kdd11);
            A12 = 2/tstep*Mdd11*(2/tstep*(ModelDx - ModelDxn_1) - ModelVxn_1) + N;
%             b1 = -(1+lener)*(4/tstep^2*Mdd11*ModelDx + N) + Mdd11*(ModelAxn_1 + 4/tstep^2*((1+lener)*ModelDxn_1 + tstep*(1+lener/2)*ModelVxn_1));
            b1 = -(Mdd11*((1+lener)*4/tstep^2*(ModelDx-ModelDxn_1) - ModelAxn_1 - 4/tstep*(1+lener/2)*ModelVxn_1) + N); % try to be better conditioned
            isw = 12; % minus sign on ModelVxn_1 term above fixed 6/17/13
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
            end
            b2 = En - SysEner;
            Fd1 = b1;
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                y = A11\b1;
                z = A11\A12;
                del_lener = (A12'*y - b2)/(A12'*z);
                del_ModelDx = y - del_lener*z;

                % Update field values with increments
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                lener = lener + del_lener;

                % Form tangent matrix and residual vector

                GetAllLoads
    %               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                    Kn = Kprime + coeffk*Kdd11;
                end

                % compute residual and energy balance equations
                N = -(Fd1 -(coeffkl2*Fext1 - F1n_1) + Mdd11*ModelAx); %flip sign, subtract off other stuff done in element routine
                A11 = (1+lener)*(4/tstep^2*Mdd11 + Kdd11);
                A12 = 2/tstep*Mdd11*(2/tstep*(ModelDx - ModelDxn_1) - ModelVxn_1) + N;
                b1 = -(1+lener)*(4/tstep^2*Mdd11*ModelDx + N) + Mdd11*(ModelAxn_1 + 4/tstep^2*((1+lener)*ModelDxn_1 + tstep*(1+lener/2)*ModelVxn_1));
                isw = 12;
                FormFE
                if numberlinear > 0
                SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
                end
                b2 = En - SysEner;
                Rnorm = norm(b1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end

        %%
        case 4
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

%             ModelDx = ModelDxn_1 + tstep*ModelVxn_1;
%             ModelVx = ModelVxn_1;
%             ModelAx = zeros(neq,1);
            ModelDx = ModelDxn_1 + coeff6*ModelVxn_1 + coeff7*ModelAxn_1;
            ModelVx = ModelVxn_1 + coeff8*ModelAxn_1;
            ModelAx = zeros(neq,1);
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Kn = Kprime + Kdd11;
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                del_ModelDx = Kn\Fd1;
                
                % Update field values with increments
                
                ModelDx = ModelDx + del_ModelDx;
                
                % Perform line search on increments
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                    Kn = Kprime + Kdd11;
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            
        %%
        case 5
            
            lamda = lamda + s_del_lamda_n1;
            
            %Set predictor values for d,v,a (d initialized as d_n-1)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                if numbernonlinear == 0 && exist('factorize','file')
%                 yvec = (KnR')\(Schol'*Fd1);
%                 xvec = KnR\yvec;
%                 del_ModelDx = Schol*xvec;
                del_ModelDx = KnR\Fd1;
                elseif numbernonlinear == 0
                del_ModelDx = Kdd11\Fd1;
                else
%                 Kn = Kprime + coeff5*Kdd11;
                del_ModelDx = Kdd11\Fd1;
                end
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                
                % Perform line search on increments
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end

        %%
        case 6
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
                     
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                del_ModelDx = Kdd11\Fd1;

                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;

                % Form tangent matrix and residual vector

                GetAllLoads
    %           Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
        
        %%
        case 7
            
            s_del_ModelDx = zeros(neq,1);
            
            %arc length parameters
            if b > 0 %force control or arc-length
            
                %Compute tangent stiffness

                isw = 3;
                FormFE%LNL
                
                detKn = det(Kdd11);
                qbar = Kdd11\(Fc1+Fc1np+FcU+Fdisp);
                f = getf(b,c,K0,qbar,1,neq);
                s_del_lamda_n1 = abs(s_del_a/f);
                %sign del_lamda
                if sign(detKn) == -sign(detKn_1) && (detKn_2 - detKn_1)*sign(detKn_1) >= 0
                    s_del_lamda_n1 = -sign(s_del_lamda_n0)*s_del_lamda_n1;
                else
                    s_del_lamda_n1 = sign(s_del_lamda_n0)*s_del_lamda_n1;
                end
                s_del_lamda_n0 = s_del_lamda_n1;
                %Increment quantities in proper direction
                s_del_ModelDx = s_del_lamda_n1*qbar; %1/17/13
                lamda = lamda + s_del_lamda_n1;
%                 linesearchFE
                s=1; %1/17/13
                ModelDx = ModelDx + s*s_del_ModelDx; %1/17/13 %%%%%% NEEDS TO BE BACK FOR NL PROBLEMS IN NL MIXED ELASTICITY
%                 COMMENTED 4/18/13 TO GET BODY FORCE LINEAR PROB TO WORK
                f = getf(b,c,K0,s_del_ModelDx,s_del_lamda_n1,neq); %1/17/13
                residtol = residratio*s_del_a;
            else %displacement control
                lamda = lamda + s_del_lamda_n1;
                residtol = s_del_a;
            end
            rarc = s_del_a;

            initia = 0;
            
            if b > -1
                f = getf(b,c,K0,s_del_ModelDx,s_del_lamda_n1,neq);
                rarc = abs(f- s_del_a);
            end
            
            GetAllLoads
            isw = 3;
            FormFE%LNL
            Kn = Kdd11;
            GetAllBCs;
            Fd1 = Fd1 - Kdf1*(gBC-gBC_n);
            
            qbar2 = Kn\Fd1;
            qbar = Kn\(Fc1+Fc1np+FcU+Fdisp);
            if b > -1 %force control or arc-length
                %Compute small delta-lamda
                dfdd = c/f*K0.*s_del_ModelDx;
                dfdl = b/f*s_del_lamda_n1;
                denomdl = dfdd'*qbar;
                denomdl = denomdl + dfdl;
                numerdl = dfdd'*qbar2;
                numerdl = s_del_a - f - numerdl;

                del_lamda = numerdl/denomdl;
                del_ModelDx = qbar2 + del_lamda*qbar;
            else %displacement control without nodal forces
                del_ModelDx = qbar2;
            end
                
            % Update field values with increments
                
            if b > -1
                lamda = lamda + del_lamda;
                if lamda > lamdamax
                    b = 1;
                    c = 0;
                    lamda = lamdamax;
                    s_del_lamda_n1 = 0;
                    stepmax = step + 1;
                    break
                end

                s_del_lamda_n1 = s_del_lamda_n1 + del_lamda;
                    
                s = 1;
                    
                del_ModelDx = s*del_ModelDx;

            end
           
            s = 1;
            s_del_ModelDx = s_del_ModelDx + s*del_ModelDx;
            ModelDx = ModelDx + s*del_ModelDx;

            if b > -1
                f = getf(b,c,K0,s_del_ModelDx,s_del_lamda_n1,neq);
                rarc = abs(f- s_del_a);
            end
                
            % Form tangent matrix and residual vector

            GetAllLoads
            GetAllBCs;
%               Kn = Kdd11;
            if mnr == 1
                isw = 6;
                FormFE

            else
                isw = 3;
                FormFE%LNL
                Kn = Kdd11;
%                     detKn = det(Kdd11);
            end
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e | %+1.7e | %+1.7e\n',Rnorm,log(Rnorm)/log(10),0)
                Rnorm1 = Rnorm;
            end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            %compute initial tangent and residual

        %%
        case 8
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

            % Predict accel as zero
            ModelDx(ddofs) = ModelDxn_1(ddofs) + coeff6*ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs);
            ModelVx(ddofs) = ModelVxn_1(ddofs) + coeff8*ModelAxn_1(ddofs);
            ModelAx = zeros(neq,1);
%             % Predict displ as zero
%             ModelDx(ddofs) = 0;
%             ModelVx(ddofs) = Ngamma*tstep*-1/(Nbeta*tstep^2)*(ModelDxn_1(ddofs) + coeff6*ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs)) + (ModelVxn_1(ddofs) + coeff8*ModelAxn_1(ddofs));
%             ModelAx(ddofs) = -1/(Nbeta*tstep^2)*(ModelDxn_1(ddofs) + coeff6*ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs));
%             % Predict veloc as zero
%             ModelDx(ddofs) = ModelDxn_1(ddofs) + coeff6*ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs) + Nbeta*tstep^2*-1/(Ngamma*tstep)*(ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs));
%             ModelVx(ddofs) = 0;
%             ModelAx(ddofs) = -1/(Ngamma*tstep)*(ModelVxn_1(ddofs) + coeff7*ModelAxn_1(ddofs));
            
%             Predict dthet as zero
            ModelDx(tdofs) = ModelDxn_1(tdofs) + (1-Nalphat)*tstep*ModelVxn_1(tdofs);
            ModelVx(tdofs) = zeros(length(tdofs),1);
%             % Predict theta as zero
%             ModelDx(tdofs) = 0;
%             ModelVx(tdofs) = -1/(Nalphat*tstep)*(ModelDxn_1(tdofs) + (1-Nalphat)*tstep*ModelVxn_1(tdofs));
                     
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                del_ModelDx = Kdd11\Fd1;

                % Update field values with increments

                ModelDx(ddofs) = ModelDx(ddofs) + coeff3*del_ModelDx(ddofs);
                ModelVx(ddofs) = ModelVx(ddofs) + coeff2*del_ModelDx(ddofs);
                ModelAx(ddofs) = ModelAx(ddofs) + coeff0*del_ModelDx(ddofs);

                ModelDx(tdofs) = ModelDx(tdofs) + del_ModelDx(tdofs);
                ModelVx(tdofs) = ModelVx(tdofs) + 1/(Nalphat*tstep)*del_ModelDx(tdofs);

                % Form tangent matrix and residual vector

                GetAllLoads
    %           Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end

        %%
        case 10 % Explicit dynamics
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 6;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
            del_ModelDx = Mdd11\Fd1;

            % Update field values with increments

            ModelDx = ModelDx + coeffd*del_ModelDx;
            ModelVx = ModelVx + coeffv*del_ModelDx;
            ModelAx = ModelAx + coeffa*del_ModelDx;
            
            Rnorm = 0;
            
        %%
        case 11
            
            s_del_ModelDx = zeros(neq,1);
            if numincrem1 > 1 % Interpolate between steps for current adapted increment
              if step > 1
                lamda = mults(step-1)+(mults(step)-mults(step-1))*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
              else
                lamda = mults(step)*((increment1-1)/numincrem1+increment2/numincrem1/numincrem2);
              end
            else
              lamda = mults(step);
            end
            residtol = 1;
            rarc = 0;

            initia = 1;
            subcycle = 1;
            maxsubcyc = 0;
            
            if step == 1 && increment1 == 1
                GetAllLoads
                isw = 3;
                FormFE
                Kn = Kdd11;
                GetAllBCs;
                Fd1 = Fd1 - Kdf1*(gBC-gBC_n);
            else
                % Extrapolation of displacements like Warp3D
                GetAllLoads
                GetAllBCs;
                ModelDx = ModelDx + ex_del_ModelDx;
                isw = 3;
                FormFE
                Kn = Kdd11;
            end
            
            initia = 0;

%             % Perform smoothing/interpolation to compute Fe_inv
%             isw = 98;
%             FormNye
            
            del_ModelDx = Kn\Fd1;
           
            s = 1;
            ModelDx = ModelDx + s*del_ModelDx;
                
            % Form tangent matrix and residual vector

            GetAllLoads
            GetAllBCs;
%               Kn = Kdd11;
            if mnr == 1
                isw = 6;
                FormFE
                if failFEelem % an element failed
                    incrementfailed = 1;
                    break % out of initial step
                end

            else
                if LieHist == 1
                    % Subcycle to update all history variables implicitly
                    subtable = 1;
                    %% continue 110
                    for subcycle = 1:numsubcycles+1
                        SubNorms = zeros(numSubnormvals,1);
                        isw = 3;
                        FormFE
                        if failFEelem % an element failed
                            incrementfailed = 1;
                            break % out of initial step
                        end
                        % Perform smoothing/interpolation to compute Fe_inv
                        isw = 100;
                        FormNyeSub
                        if subcycle > 1
                            SubNorms(1) = norm(reshape(LieList_old-LieList,numLie*numnp*nummat,1));
                            SubnormsL(1:3,subtable,step) = SubNorms(1:3); % absolute norms
                            SubnormsL(4:6,subtable,step) = SubNorms(1:3)./SubNorm(1:3); % relative norms
                            subtable = subtable + 1;
                            if SubNorms(1) < subcycTol*SubNorm(1) && ...
                               SubNorms(2) < subcycTol*SubNorm(2) && ...
                               SubNorms(3) < subcycTol*SubNorm(3)
                                break
                            end
                        end
                        LieList_old = LieList(1:numLie,:,:);
                        if subcycle == 1
                            SubNorm = zeros(3,1);
                            SubNorm(1) = norm(reshape(LieList(1:numLie,1:numnp,1:nummat),numLie*numnp*nummat,1));
                            SubNorm(2) = SubNorms(2); % stress
                            SubNorm(3) = SubNorms(3); % tt
                            SubnormsT(1:3,iter+1,step) = SubNorm(1:3);
                        end
                    end % go to 110
                    if failFEelem % an element failed
                        break % out of initial step
                    end
                    if printRnorm
                    subcycle
                    end
                    maxsubcyc = max(maxsubcyc,subcycle);
                else
                    isw = 3;
                    FormFE
                    if failFEelem % an element failed
                        incrementfailed = 1;
                        break % out of initial step
                    end
                end
                Kn = Kdd11;
%                     detKn = det(Kdd11);
            end

%             % Perform smoothing/interpolation to compute Fe_inv
%             isw = 98;
%             FormNye
            
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e  | %1.7e  | %+1.7e  | %1.7e\n',Rnorm,norm(Fd1,inf),log(Rnorm)/log(10),0)
                Rnorm1 = Rnorm;
            end
            Rnorm1 = Rnorm;
            Rnorm2 = Rnorm;
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            
        case 12
            
            lamda = lamda + s_del_lamda_n1;
    
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 6;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
            del_ModelDx = Mdd11\Fd1;

            % Update field values with increments

            ModelDx = ModelDx + coeffd*del_ModelDx;
            ModelVx = ModelVx + coeffv*del_ModelDx;
            ModelAx = ModelAx + coeffa*del_ModelDx;
            
            Rnorm = 0;
            
        case 13 % GLS
            
            initia = 1;
            lamda = lamda + s_del_lamda_n1;
            
            %Set predictor values for d,v,a (a initialized as zero)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            if enercons == 1 % initialize LM for energy to zero
                lener = 0;
            end
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            % first initialize assuming lener=0, meaning linear elasticity

            GetAllLoads
            isw = 3;
            FormFE
                    
            A11 = Kdd11;%(1+lener)*(4/tstep^2*Mdd11 + Kdd11);
            b1 = Fd1;
            y = A11\b1;
            z = 0*y;
            del_ModelDx = y;
            del_lener = 0;
            ModelDx = ModelDx + coeffd*del_ModelDx;
            ModelVx = ModelVx + coeffv*del_ModelDx;
            ModelAx = ModelAx + coeffa*del_ModelDx;

            initia = 0;
            %compute initial tangent and residual, lener is still 0

            GetAllLoads
            isw = 3;
            FormFE
            Kn = Kprime + coeffk*Kdd11;
            
            A11 = Kdd11;
            b1 = Fd1;
            isw = 80;
            FormA12
            A12 = ModelA12 + Fext1;
            isw = 12;
            FormFE
            if numberlinear > 0
            SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
            end
            b2 = En - SysEner;
            isw = 6;
            FormFE
            A21 = (2/tstep*Mdd11*(2/tstep*(ModelDx - ModelDxn_1) - ModelVxn_1) - Fd1)';
            Fd1 = b1;
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                y = A11\b1;
                z = A11\A12;
                del_lener = (A21*y - b2)/(A21*z);
                del_ModelDx = y - del_lener*z;

                % Update field values with increments
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                lener = lener + del_lener;

                GetAllLoads
                isw = 3;
                FormFE
                Kn = Kprime + coeffk*Kdd11;

                A11 = Kdd11;
                b1 = Fd1;
                isw = 80;
                FormA12
                A12 = ModelA12 + Fext1;
                isw = 12;
                FormFE
                if numberlinear > 0
                SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
                end
                b2 = En - SysEner;
                isw = 6;
                FormFE
                A21 = (2/tstep*Mdd11*(2/tstep*(ModelDx - ModelDxn_1) - ModelVxn_1) - Fd1)';
                Fd1 = b1;
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            
        case 14
            
            lamda = lamda + s_del_lamda_n1;
            
            %Set predictor values for d,v,a (d initialized as d_n-1)

            ModelDx = coeffdd*ModelDxn_1 + coeffdv*ModelVxn_1 + coeffda*ModelAxn_1;
            ModelVx = coeffvd*ModelDxn_1 + coeffvv*ModelVxn_1 + coeffva*ModelAxn_1;
            ModelAx = coeffad*ModelDxn_1 + coeffav*ModelVxn_1 + coeffaa*ModelAxn_1;
            
            %Set arc-length parameters so as not to affect time integration

            s_del_a = 1;
            residtol = s_del_a;
            rarc = 0.9*s_del_a;

            %compute initial tangent and residual

            GetAllLoads
            isw = 3;
            FormFE
            Rnorm = norm(Fd1);
            if printRnorm
                fprintf('%1.7e\n',Rnorm)
            end
                
                if numbernonlinear == 0
%                 yvec = (KnR')\(Schol'*Fd1);
%                 xvec = KnR\yvec;
%                 del_ModelDx = Schol*xvec;
                del_ModelDx = KnR\Fd1;
                else
%                 Kn = Kprime + coeff5*Kdd11;
                del_ModelDx = Kdd11\Fd1;
                end
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                
                % Perform line search on increments
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
            if IHist == 1
                IterList(iter+4,step) = Rnorm;
            end
            
    end % Initialization switch, line 50
    
    %% Set residual tolerance
    if Rnorm < 50*Residratio %already converged, i.e. linear
    Residtol = 1.2*Rnorm;
    if transient == 7
    residtol = 1.2*rarc;
    end
    else
    Residtol = max(Residratio*Rnorm,40*eps);
    end
    

    
    
    %% continue 60
    % Loop for Newton-Raphson Method NL Solver
    while iter < itermax && (Rnorm > Residtol || (transient==7 && rarc > residtol)) && numbernonlinear > 0
        
        iter = iter + 1;
        
        %% switch line 70
        % Solve for incremental quantities
        switch transient
            
            %%
            case 0
                
                del_ModelDx = Kn\Fd1;
                
                % Update field values with increments
                
                % Perform line search
                linesearchFE

                ModelDx = ModelDx + sLS*del_ModelDx;
                
                % Form tangent matrix and residual vector
         
                GetAllLoads
                GetAllBCs;
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                    if failFEelem % an element failed
                        break % out of initial step
                    end
                else
                    isw = 3;
                    FormFE
                    if failFEelem % an element failed
                        break % out of initial step
                    end
                    Kn = Kdd11;
%                     detKn = det(Kdd11);
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e | %+1.7e | %+1.7e\n',Rnorm,log(Rnorm)/log(10),log(Rnorm)/log(10) - 2*log(Rnorm1)/log(10))
                    Rnorm1 = Rnorm;
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end

            %%
            case {1,2}
                
                if numbernonlinear == 0 && exist('factorize','file')
%                 yvec = (KnR')\(Schol'*Fd1);
%                 xvec = KnR\yvec;
%                 del_ModelDx = Schol*xvec;
                del_ModelDx = KnR\Fd1;
                elseif numbernonlinear == 0
                del_ModelDx = Kdd11\Fd1;
                else
%                 Kn = Kprime + coeff5*Kdd11;
                del_ModelDx = Kdd11\Fd1;
                end
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                
                % Perform line search on increments
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end
                
            %%
            case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
                
                y = A11\b1;
                z = A11\A12;
                del_lener = (A12'*y - b2)/(A12'*z);
                del_ModelDx = y - del_lener*z;
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
                lener = lener + del_lener;
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                    Kn = Kprime + coeffk*Kdd11;
                end
                    
                % compute residual and energy balance equations
                N = -(Fd1 -(coeffkl2*Fext1 - F1n_1) + Mdd11*ModelAx); %flip sign, subtract off other stuff done in element routine
                A11 = (1+lener)*(4/tstep^2*Mdd11 + Kdd11);
                A12 = 2/tstep*Mdd11*(2/tstep*(ModelDx - ModelDxn_1) - ModelVxn_1) + N;
                b1 = -(1+lener)*(4/tstep^2*Mdd11*ModelDx + N) + Mdd11*(ModelAxn_1 + 4/tstep^2*((1+lener)*ModelDxn_1 + tstep*(1+lener/2)*ModelVxn_1));
                isw = 12;
                FormFE
                if numberlinear > 0
                SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
                end
                b2 = En - SysEner;
                Rnorm = norm(b1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end

            %%
            case 4
                
                del_ModelDx = Kn\Fd1;
                
                % Update field values with increments
                
                ModelDx = ModelDx + del_ModelDx;
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                    Kn = Kprime + Kdd11;
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end

            %%
            case 5
                
                if numbernonlinear == 0 && exist('factorize','file')
%                 yvec = (KnR')\(Schol'*Fd1);
%                 xvec = KnR\yvec;
%                 del_ModelDx = Schol*xvec;
                del_ModelDx = KnR\Fd1;
                elseif numbernonlinear == 0
                del_ModelDx = Kdd11\Fd1;
                else
%                 Kn = Kprime + coeff5*Kdd11;
                del_ModelDx = Kdd11\Fd1;
                end
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end

            %%
            case 6
                
                del_ModelDx = Kdd11\Fd1;
                
                % Update field values with increments
                
                ModelDx = ModelDx + coeffd*del_ModelDx;
                ModelVx = ModelVx + coeffv*del_ModelDx;
                ModelAx = ModelAx + coeffa*del_ModelDx;
        
                % Form tangent matrix and residual vector
         
                GetAllLoads
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end
            
            %%
            case 7
                
                qbar2 = Kn\Fd1;
                qbar = Kn\(Fc1+Fc1np+FcU+Fdisp);
                if b > -1 %force control or arc-length
                    %Compute small delta-lamda
                    dfdd = c/f*K0.*s_del_ModelDx;
                    dfdl = b/f*s_del_lamda_n1;
                    denomdl = dfdd'*qbar;
                    denomdl = denomdl + dfdl;
                    numerdl = dfdd'*qbar2;
                    numerdl = s_del_a - f - numerdl;

                    del_lamda = numerdl/denomdl;
                    del_ModelDx = qbar2 + del_lamda*qbar;
                else %displacement control without nodal forces
                    del_ModelDx = qbar2;
                end
                
                % Update field values with increments
                
                if b > -1
                    lamda = lamda + del_lamda;
                    if lamda > lamdamax
                        b = 1;
                        c = 0;
                        lamda = lamdamax;
                        s_del_lamda_n1 = 0;
                        stepmax = step + 1;
                        break
                    end

                    s_del_lamda_n1 = s_del_lamda_n1 + del_lamda;
                    
                % Perform line search on increments
% %                     if iter >= 2
% %                         smax = 1;
% %                     end
% %                     s = smax;
% if lamda >= 2
%     if lamda == 2.5
% if iter < 20
% %                     linesearchFE
%                     s = 0.2;
% else
%     s = 1;
% end
%     else
% if iter < 3
% %                     linesearchFE
%                     s = 0.5;
% else
    s = 1;
% end
%     end
% else
%     s = 1;
% end
                    
                    del_ModelDx = s*del_ModelDx;

                end
                
% if numdbond > 0
%     linesearchFE2
% %     s
% else
    s = 1;
% end
                s_del_ModelDx = s_del_ModelDx + s*del_ModelDx;
                ModelDx = ModelDx + s*del_ModelDx;
                
                if b > -1
                    f = getf(b,c,K0,s_del_ModelDx,s_del_lamda_n1,neq);
                    rarc = abs(f- s_del_a);
                end
                
                % Form tangent matrix and residual vector
         
                GetAllLoads
                GetAllBCs;
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                    
                else
                    isw = 3;
                    FormFE
                    Kn = Kdd11;
%                     detKn = det(Kdd11);
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e | %+1.7e | %+1.7e\n',Rnorm,log(Rnorm)/log(10),log(Rnorm)/log(10) - 2*log(Rnorm1)/log(10))
                    Rnorm1 = Rnorm;
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end
                
            %%
            case 8
                
                del_ModelDx = Kdd11\Fd1;

                % Update field values with increments

                ModelDx(ddofs) = ModelDx(ddofs) + coeff3*del_ModelDx(ddofs);
                ModelVx(ddofs) = ModelVx(ddofs) + coeff2*del_ModelDx(ddofs);
                ModelAx(ddofs) = ModelAx(ddofs) + coeff0*del_ModelDx(ddofs);

                ModelDx(tdofs) = ModelDx(tdofs) + del_ModelDx(tdofs);
                ModelVx(tdofs) = ModelVx(tdofs) + 1/(Nalphat*tstep)*del_ModelDx(tdofs);

                % Form tangent matrix and residual vector

                GetAllLoads
    %           Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                else
                    isw = 3;
                    FormFE
                end
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e\n',Rnorm)
                end
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end
            
            %%
            case 11
                
                del_ModelDx = Kn\Fd1;
                
                % Update field values with increments
                
% if numdbond > 0
%     linesearchFE2
% %     s
% else
    s = 1;
% end
                ModelDx = ModelDx + s*del_ModelDx;
                
                % Form tangent matrix and residual vector
         
                GetAllLoads
                GetAllBCs;
%               Kn = Kdd11;
                if mnr == 1
                    isw = 6;
                    FormFE
                    if failFEelem % an element failed
                        break % out of NR loop
                    end
                    
                else
                    if LieHist == 1
                        %% continue 120
                        % Subcycle to update all history variables implicitly
                        for subcycle = 1:numsubcycles+1
                            SubNorms = zeros(numSubnormvals,1);
                            isw = 3;
                            FormFE
                            if failFEelem % an element failed
                                break % out of NR loop
                            end
                            % Perform smoothing/interpolation to compute Fe_inv
                            isw = 100;
                            FormNyeSub
                            if subcycle > 1
                                SubNorms(1) = norm(reshape(LieList_old-LieList,numLie*numnp*nummat,1));
                                SubnormsL(1:3,subtable,step) = SubNorms(1:3); % absolute norms
                                SubnormsL(4:6,subtable,step) = SubNorms(1:3)./SubNorm(1:3); % relative norms
                                subtable = subtable + 1;
                                if SubNorms(1) < subcycTol*SubNorm(1) && ...
                                   SubNorms(2) < subcycTol*SubNorm(2) && ...
                                   SubNorms(3) < subcycTol*SubNorm(3)
                                    break
                                end
                            end
                            LieList_old = LieList(1:numLie,:,:);
                            if subcycle == 1
                                SubNorm = zeros(3,1);
                                SubNorm(1) = norm(reshape(LieList,numLie*numnp*nummat,1));
                                SubNorm(2) = SubNorms(2); % stress
                                SubNorm(3) = SubNorms(3); % tt
                                SubnormsT(1:3,iter+1,step) = SubNorm(1:3);
                            end
                        end % go to 120
                        if failFEelem % an element failed
                            break % out of NR loop
                        end
                        if printRnorm
                        subcycle
                        end
                        maxsubcyc = max(maxsubcyc,subcycle);
                    else
                        isw = 3;
                        FormFE
                        if failFEelem % an element failed
                            break % out of NR loop
                        end
                    end
                    Kn = Kdd11;
%                     detKn = det(Kdd11);
                end

%                 % Perform smoothing/interpolation to compute Fe_inv
%                 isw = 98;
%                 FormNye
                
                Rnorm = norm(Fd1);
                if printRnorm
                    fprintf('%1.7e  | %1.7e  | %+1.7e  | %1.7e\n',Rnorm,norm(Fd1,inf),log(Rnorm)/log(10),log(Rnorm)/log(10) - 2*log(Rnorm1)/log(10))
                end
                if Rnorm > 5*Rnorm1 && Rnorm > 5*Rnorm2;
                    diverge_flag = 1;
                    break
                end
                Rnorm2 = Rnorm1;
                Rnorm1 = Rnorm;
                if IHist == 1
                    IterList(iter+4,step) = Rnorm;
                end

        end
        
    end % go to 60, End NRM NL Solver
    
    
    
    
    %% Check for convergence; code comes to this line if an element failed
    % in an iteration of NR
    incrementfailed = ((Rnorm > Residtol || (transient==7 && rarc > residtol)) && numbernonlinear > 0); % successful convergence
    incrementfailed = incrementfailed && Rnorm > 1e-3; % catch iterations stuck at machine epsilon
    incrementfailed = incrementfailed || failFEelem; % also find if an element failed
    incrementfailed = incrementfailed || diverge_flag; % also find if simulation diverging
    
    
    %% Update the solution fields for an increment, only if not last increment
    if ~incrementfailed && (increment1 < numincrem1 || increment2 < numincrem2)
    
        
      %% Perform data shifts that are global and should be done prior to post-process
      switch transient
            
        case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
            
            ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
            ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
            
        case 4
            
            ModelVx = 2/tstep*(ModelDx-ModelDxn_1) - ModelVxn_1;
            ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
            
%         case 13 % GLS
%             
%             ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
%             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
            
    end
    
      %% switch 70
      % Perform data shifts required after post-process
      switch transient
            
        case 0
            
            %increment values for step n+1
            ModelDxn_1 = ModelDx;
%             s_del_lamda_n0 = s_del_lamda_n1;
            gBC_n = gBC;
            
        case {1,2}
            
 %            ModelDx = (ModelDx - ModelDxn_1)/(1+alpha) + ModelDxn_1;
 %            ModelVx = (ModelVx - ModelVxn_1)/(1+alpha) + ModelVxn_1;
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            
        case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
            
%             ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
%             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            gBC_n = gBC;
            
        case 4
            
%             ModelVx = 2/tstep*(ModelDx-ModelDxn_1) - ModelVxn_1;
%             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            gBC_n = gBC;
            
        case 5
            
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            gBC_n = gBC;
            
        case 6
            
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            gBC_n = gBC;
            
        case 7
            
            %increment values for step n+1
            if b > -1
                detKn_2 = detKn_1;
                detKn_1 = detKn;
            end
            ModelDxn_1 = ModelDx;
%             s_del_lamda_n0 = s_del_lamda_n1;
            gBC_n = gBC;
            
        case 8
            
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;
            gBC_n = gBC;
                    
        case 10
            
 %            ModelDx = (ModelDx - ModelDxn_1)/(1+alpha) + ModelDxn_1;
 %            ModelVx = (ModelVx - ModelVxn_1)/(1+alpha) + ModelVxn_1;
            ModelDxn_1 = ModelDx;
            ModelVxn_1 = ModelVx;
            ModelAxn_1 = ModelAx;    
            
        case 11
            
            %increment values for step n+1
            ex_del_ModelDx = ModelDx - ModelDxn_1; % extrapolation displacement
            ModelDxn_1 = ModelDx;
%             s_del_lamda_n0 = s_del_lamda_n1;
            gBC_n = gBC;
            
%             if LieHist == 1
%                 % Perform smoothing/interpolation to compute Fe_inv
%                 isw = 98;
%                 FormNye
%             end
      end % Data shifts, line 70
        
      %% Update internal variables
      hrn1 = 2;
      hrn2 = 1;
      reshis
        
    else % increment failed or final increment; exit and subdivide
        
      break
        
    end % if, Update the solution fields for an increment
    
    
   end % go to 40, increment2
    
    
   %% Reset time step from level 2 increment
   if incrementfailed
       break % skip and exit out to main checks
   elseif ~(increment1 < numincrem1 || increment2 < numincrem2)
       % multiply time step
       tstep = tstep*numincrem2;
       if exist('ex_del_ModelDx','var')
           ex_del_ModelDx = ex_del_ModelDx*numincrem2;
       end
       % reset number of increments on level 2
       numincrem2 = 1;
   end
    
    
  end % go to 30, increment1
   
   
  %% NR loop has failed or failFEelem==1; so, adapt the time step
  % Code comes to this point if element fails during initialization
  incrementfailed = incrementfailed || failFEelem; % also find if an element failed
    
  if incrementfailed && numincrem1 == 1 && numincrem2 == 1 % step failed, subdivide once
        
      stepnotdone1 = 1;
      disp('subdividing first time')
      numincrem1 = 4;
      % divide time step
      tstep = tstep/numincrem1;
      % reset solution to last converged increment/step
      ModelDx = ModelDxn_1;
      ModelVx = ModelVxn_1;
      ModelAx = ModelAxn_1;
      if exist('ex_del_ModelDx','var')
          ex_del_ModelDx = ex_del_ModelDx/numincrem1;
      end
      %% Reset internal variables
      hrn1 = 1;
      hrn2 = 2;
      temphrstore = hrstore;
      hrstore = 0;
      reshis
      hrstore = temphrstore;
        
  elseif incrementfailed && numincrem1 == 4 && numincrem2 == 1 % increment failed, subdivide twice
      
      stepnotdone1 = 1;
      disp('subdividing second time')
      incremstart1 = increment1;
      numincrem2 = 4;
      % divide time step
      tstep = tstep/numincrem2;
      % reset solution to last converged increment/step
      ModelDx = ModelDxn_1;
      ModelVx = ModelVxn_1;
      ModelAx = ModelAxn_1;
      if exist('ex_del_ModelDx','var')
          ex_del_ModelDx = ex_del_ModelDx/numincrem2;
      end
      %% Reset internal variables
      hrn1 = 1;
      hrn2 = 2;
      temphrstore = hrstore;
      hrstore = 0;
      reshis
      hrstore = temphrstore;
        
  elseif incrementfailed && numincrem1 == 4 && numincrem2 == 4 % terminate the simulation
      
      error('Failed to converge after subdividing twice')
      
  else % successful step, end and move on
      
      stepnotdone1 = 0;
      tstep = tstep_step;
      %         ex_del_ModelDx = ex_del_ModelDx*numincrem2; % handled at end instead
      
  end
    
    
  end % go to 20, stepnotdone
    
    
  
  
 %% Perform data shifts that are global and should be done prior to post-process
 switch transient
      
     case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
          
         ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
         ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
          
     case 4
          
         ModelVx = 2/tstep*(ModelDx-ModelDxn_1) - ModelVxn_1;
         ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
          
 end
    
    
 %% Output converged values to history arrays
 TimeList(step) = lamda;
 if IHist == 1
     IterList(1,step) = iter;
     if iter < itermax
         IterList(2,step) = 1;
     else
         IterList(2,step) = Rnorm/Residtol;
     end
     if LieHist == 1
         IterList(3,step) = maxsubcyc;
     end
 end
 %% Displacements
 if DHist == 1
     DispList(DOFa+(step-1)*ndf*numnp) = ModelDx(NDOFT2(DOFa));
     DispList(DOFi+(step-1)*ndf*numnp) = gBC(NDOFT2(DOFi)-neq);
 end
 %% Velocities
 if VHist == 1
     VeloList(DOFa+(step-1)*ndf*numnp) = ModelVx(NDOFT2(DOFa));
 end
 %% Accelerations
 if AHist == 1
     AcceList(DOFa+(step-1)*ndf*numnp) = ModelAx(NDOFT2(DOFa));
 end
 %% Reaction forces
 if FHist == 1
     Fext1 = zeros(neq,1);
     if LieHist == 1
         %% continue 130
         % Subcycle to update all history variables implicitly
         for subcycle = 1:numsubcycles+1
             SubNorms = zeros(numSubnormvals,1);
             isw = 6;
             FormFE
             % Perform smoothing/interpolation to compute Fe_inv
             isw = 100;
             FormNyeSub
             if subcycle > 1
                 SubNorms(1) = norm(reshape(LieList_old-LieList,numLie*numnp*nummat,1));
                 if SubNorms(1) < subcycTol*SubNorm(1) && ...
                         SubNorms(2) < subcycTol*SubNorm(2) && ...
                         SubNorms(3) < subcycTol*SubNorm(3)
                     break
                 end
             end
             LieList_old = LieList(1:numLie,:,:);
             if subcycle == 1
                 SubNorm = zeros(3,1);
                 SubNorm(1) = norm(reshape(LieList,numLie*numnp*nummat,1));
                 SubNorm(2) = SubNorms(2); % stress
                 SubNorm(3) = SubNorms(3); % tt
             end
         end % go to 130
     else
         isw = 6;
         FormFE
     end
     
     ForcList(DOFa+(step-1)*ndf*numnp) = Fd1(NDOFT2(DOFa));
     ForcList(DOFi+(step-1)*ndf*numnp) = Fd3(NDOFT2(DOFi)-neq);
 end
 
 %% Process nodal stresses, only at desired step in SHList
 if SHist == 1 && stepS <= length(SHList) && step == SHList(stepS)
     
     StreList2 = zeros(numnp,npstr);
     Eareas = zeros(numnp,1);
     
     isw = 25;
     FormS2
     
     for stres = 1:npstr
         StreList2(:,stres) = StreList2(:,stres)./Eareas;
     end
     StreList(1:npstr,:,stepS) = StreList2';
     
     stepS = stepS + 1;
     
 end
 
 %% Process elemental stresses, only at desired step in SHList
 if SEHist == 1 && stepSE <= length(SEHList) && step == SEHList(stepSE)
     
     StreList2 = zeros(numel,nestr);
     
     isw = 26;
     FormFE
     
     StreListE(1:nestr,:,stepSE) = StreList2';
     
     stepSE = stepSE + 1;
     
 end
 
 %% Plasticity variables
 if PHist == 1
     isw = 24;
     FormFE
 end
 %% Macroscopic stress/strain for RVE
 if SSHist == 1
     isw = 51;
     FormFE
     
     SSList(1:12,step) = SSValues(1:12)/SSValues(13);
     SSList(13,step) = SSValues(13);
 end
 %% Interface segment quantities/fields
 if IFHist == 1
     isw = 60;
     FormI
 end
 
 %% J-Integrals
 if JHist == 1
     isw = 16;
     JIntegral
     JintList(1:ndm,1:JDImax,step) = JstepList;
 end
 
 %% Energy
 switch transient
     
     case 0
         
         if PDHist == 1
             isw = 13;
             FormFE
             PlasDiss
             PDList(step) = PlasDiss;
         end
         
     case {1,2}
         
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
                 SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
             end
             SysEner
             %                 En = SysEner;
         end
         
     case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
         
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
                 SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
             end
             SysEner
             En = SysEner
         end
         
     case 4
         
         if enercons == 1
             isw = 12;
             FormFE
             SysEner = SysEner + 1/2*ModelVx'*Mdd11*ModelVx
             %                 En = SysEner;
         end
         
     case 5
         
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
                 SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
             end
             SysEner
             %                 En = SysEner;
         end
         
     case 6
         
         if enercons == 1
             isw = 12;
             FormFE
             % For external loads; assumes trapezoidal rule
             Fext1 = Fc1np + lamda*Fc1;
             lamda_n = lamda - s_del_lamda_n1;
             Fext1_n = Fc1np + lamda_n*Fc1;
             SysEner = SysEner - 1/2*ModelDx'*(Fext1 + Fext1_n);
             SysEner
             maxen = max(maxen,SysEner);
             %                 En = SysEner;
         end
         
     case 7
         
     case 8
         
         if enercons == 1
             isw = 12;
             FormFE
             SysEner
             maxen = max(maxen,SysEner);
             %                 En = SysEner;
         end
         
     case 10
         
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
                 SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
             end
             SysEner
             %                 En = SysEner;
         end
         
     case 11
         
         %                 % Perform smoothing/interpolation to compute Fe_inv
         %                 isw = 98;
         %                 FormNye
         
                    
     case 12
        
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
             SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
             end
             SysEner
%              En = SysEner;
         end 
            
     case 13 % GLS
        
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
             SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx;
             end
             En = SysEner
         end 
            
     case 14
        
         if enercons == 1
             isw = 12;
             FormFE
             if numberlinear > 0
             SysEner = SysEner + 1/2*ModelVx'*MddLL*ModelVx + 1/2*ModelDx'*KddLL*ModelDx;
             end
             SysEner
%              En = SysEner;
         end

 end
 
 
 %% switch 80
 % Perform data shifts required after post-process
 switch transient
     
     case 0
         
         %increment values for step n+1
         ModelDxn_1 = ModelDx;
         %             s_del_lamda_n0 = s_del_lamda_n1;
         gBC_n = gBC;
         
     case {1,2}
         
         %            ModelDx = (ModelDx - ModelDxn_1)/(1+alpha) + ModelDxn_1;
         %            ModelVx = (ModelVx - ModelVxn_1)/(1+alpha) + ModelVxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         
     case 3 %Energy Conservation by LM according to Hughes, Trans of ASME
         
         %             ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
         %             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
         
     case 4
         
         %             ModelVx = 2/tstep*(ModelDx-ModelDxn_1) - ModelVxn_1;
         %             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
         
     case 5
         
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
         
     case 6
         
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
         
     case 7
         
         %increment values for step n+1
         if b > -1
             detKn_2 = detKn_1;
             detKn_1 = detKn;
         end
         ModelDxn_1 = ModelDx;
         %             s_del_lamda_n0 = s_del_lamda_n1;
         gBC_n = gBC;
         
     case 8
         
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
         
     case 10
         
         %            ModelDx = (ModelDx - ModelDxn_1)/(1+alpha) + ModelDxn_1;
         %            ModelVx = (ModelVx - ModelVxn_1)/(1+alpha) + ModelVxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         
     case 11
         
         %increment values for step n+1
         %             ex_del_ModelDx = ModelDx - ModelDxn_1; % extrapolation displacement
         ex_del_ModelDx = ModelDx - ModelDxn_1_step; % extrapolation displacement
         ModelDxn_1 = ModelDx;
         %             s_del_lamda_n0 = s_del_lamda_n1;
         gBC_n = gBC;
         
         %             if LieHist == 1
         %                 % Perform smoothing/interpolation to compute Fe_inv
         %                 isw = 98;
         %                 FormNye
         %             end

     case 12
            
 %            ModelDx = (ModelDx - ModelDxn_1)/(1+alpha) + ModelDxn_1;
 %            ModelVx = (ModelVx - ModelVxn_1)/(1+alpha) + ModelVxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx; 
            
     case 13 % GLS
            
%             ModelVx = 2/tstep*ModelDx - 2/tstep*ModelDxn_1 - ModelVxn_1;
%             ModelAx = 2/tstep*ModelVx - 2/tstep*ModelVxn_1 - ModelAxn_1;
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
            
     case 14
            
         ModelDxn_1 = ModelDx;
         ModelVxn_1 = ModelVx;
         ModelAxn_1 = ModelAx;
         gBC_n = gBC;
 end
 
 %% Store energy
 if EHist == 1 && ismember(transient,transalgos)
     EnerList(step) = SysEner;
 end
 if PDHist == 1
     PDList(step) = PlasDiss;
 end
 
 %% Perform error estimation
 if expliciterr == 1 && stepEE <= length(EEHList) && step == EEHList(stepEE)
     
     Explicit_Error_Estimation
     EEList(:,stepEE) = Energy;
     IeffList(:,:,stepEE) = Ieffvals';
     
     stepEE = stepEE + 1;
     
 end
 if impliciterr == 1 && stepIE <= length(IEHList) && step == IEHList(stepIE)
     
     Implicit_Error_Estimation
     
 end
 
 
 %% Update internal variables
 hrn1 = 2;
 hrn2 = 1;
 reshis
 
 
 %% Dump data to restart file
 if mod(step,reststep) == 0 && restartmat == 1
     if step < 10
         rfile = strcat('Restart000',num2str(step),'.mat');
     elseif step < 100
         rfile = strcat('Restart00',num2str(step),'.mat');
     elseif step < 1000
         rfile = strcat('Restart0',num2str(step),'.mat');
     else
         rfile = strcat('Restart',num2str(step),'.mat');
     end
     save(rfile)
     disp('Restart file saved!')
 end
 
 
 %% Final catch of failed load step
 if Rnorm > 1e5*Residtol
     error('Rnorm > 1 : global iterations failed to reasonably converge')
 end

end % go to 10, End time-history solution
