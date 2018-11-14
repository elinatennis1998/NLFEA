% Tim Truster
% 01/16/2014
%
% Stress-strain plot in shear for double-slip crystal plasticity example
% from deSouza text, ch 16 p. 21
% Employs the model from HYPLAS, recoded in Matlab

% values are close but don't match deSouza exactly; DEXPMP is not the
% issue, I compared that already.

clear
% clc

G = 21.10;
K = 49.98;
% tstep = 100;1000;10;
eta_max = 3;
% gmax = 10;
% gams = (0:gmax/tstep:gmax)';
% taus = 0.06+0.048*(1-exp(-gams/0.929))+0.001*gams;
% Hlist = reshape([gams'; taus'],1,2*(tstep+1));
% IPROPS = [0 0 tstep+1];
Hlist = [  0.0   0.06
  0.02  0.0693170526
  0.04116  0.0772217724
  0.06354728  0.0838439672
  0.0872330222  0.0893183217
  0.112292538  0.093780868
  0.138805505  0.0973655998
  0.166856224  0.10020138
  0.196533885  0.102409276
  0.22793285  0.104100433
  0.261152956  0.105374545
  0.296299827  0.106318961
  0.333485217  0.107008398
  0.37282736  0.107505215
  0.414451346  0.107860159
  0.458489525  0.108113452
  0.505081917  0.108296126
  0.554376668  0.108431461
  0.606530515  0.108536418
  0.661709285  0.108622998
  0.720088423  0.108699438
  0.847201058  0.108841945
  0.989486365  0.10898835
  1.14875542  0.109148551
  1.32703546  0.109327005
  1.74997629  0.109749976
  2.2799109  0.110279911
  2.94390612  0.110943906
  3.77587605  0.111775876
  5.11777643  0.113117776]';
Hlist = reshape(Hlist,1,60);
IPROPS = [0 0 30];
theta = pi/3; beta = pi/3;
NTYPE = 2;
RPROPS = [0 G K theta beta Hlist];

nstep1 = 10;
inc1 = eta_max/nstep1;
RSTAVA = [1 0 0 1 0]';
FINCR = [1 inc1; 0 1];

Stresses1 = zeros(4,nstep1);
Strains1 = zeros(5,nstep1);
for i = 1:nstep1
[DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(FINCR,IPROPS,NTYPE,RPROPS,RSTAVA);
Stresses1(:,i) = STRES;
Strains1(:,i) = RSTAVA;
end

nstep2 = 30;
inc2 = eta_max/nstep2;
RSTAVA = [1 0 0 1 0]';
FINCR = [1 inc2; 0 1];

Stresses2 = zeros(4,nstep2);
Strains2 = zeros(5,nstep2);
for i = 1:nstep2
[DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(FINCR,IPROPS,NTYPE,RPROPS,RSTAVA);
Stresses2(:,i) = STRES;
Strains2(:,i) = RSTAVA;
end

nstep3 = 300;
inc3 = eta_max/nstep3;
RSTAVA = [1 0 0 1 0]';
FINCR = [1 inc3; 0 1];

Stresses3 = zeros(4,nstep3);
Strains3 = zeros(5,nstep3);
for i = 1:nstep3
    RSTAVN = RSTAVA;
[DGAM,LALGVA,RSTAVA,STRES] = SUPDSC(FINCR,IPROPS,NTYPE,RPROPS,RSTAVA);
Stresses3(:,i) = STRES;
Strains3(:,i) = RSTAVA;
end
AMATX = CSTPDS(DGAM,1,FINCR,IPROPS,LALGVA,NTYPE,RPROPS,RSTAVA,RSTAVN,STRES);

figure(2)
clf(2)
subplot(2,2,1)
% hardening variable
hold on
plot(inc1:inc1:inc1*nstep1,Strains1(5,:),'d')
plot(inc2:inc2:inc2*nstep2,Strains2(5,:),'r+')
plot(inc3:inc3:inc3*nstep3,Strains3(5,:),'k')
title('hardening variable')
% txx
subplot(2,2,2)
hold on
plot(inc1:inc1:inc1*nstep1,Stresses1(1,:),'d')
plot(inc2:inc2:inc2*nstep2,Stresses2(1,:),'r+')
plot(inc3:inc3:inc3*nstep3,Stresses3(1,:),'k')
title('txx')
% tyy
subplot(2,2,3)
hold on
plot(inc1:inc1:inc1*nstep1,Stresses1(2,:),'d')
plot(inc2:inc2:inc2*nstep2,Stresses2(2,:),'r+')
plot(inc3:inc3:inc3*nstep3,Stresses3(2,:),'k')
title('tyy')
% txy
subplot(2,2,4)
hold on
plot(inc1:inc1:inc1*nstep1,Stresses1(3,:),'d')
plot(inc2:inc2:inc2*nstep2,Stresses2(3,:),'r+')
plot(inc3:inc3:inc3*nstep3,Stresses3(3,:),'k')
title('txy')