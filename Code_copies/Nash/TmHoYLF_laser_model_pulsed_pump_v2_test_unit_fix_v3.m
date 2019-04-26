% Based on paper "Spectroscopy and Modeling of Solid State Lanthanide
% Lasers..." (Walsh, Barnes, Petros, Yu, Singh)

% pulsed_pump_v1: Defines pump scalar pump rate, uses conditional statement in rate eq
% m-file to determine if t<=tp, in which case the pump rate is set to
% defined pump rate or 0

% pulsed_pump_v2: Defines pump pulse explicitly as a function of time, rate eq m-file
% interpolates to determine pump rate at each particular time step

% The two methods give identical results, v2 takes approx 2x longer to run

% test_unit_fix:  Fixes unit inconsistency in equations from reference
% paper (adds c*l/Lopt factor to stimulated emission terms of n7 and n8
% rate eqs, drops same factor from spontaneous emission term of photon
% density

% test_unit_fix_v2:  Same as above but uses factor of c/n instead of
% c*l/Lopt as in reference "Computational model for operation of 2um..."
% (Louchev, Urata, Saito, Wada, this paper also references original
% refernce paper "Spectroscopy and Modeling of Solid State Lanthanide
% Lasers...")

% test_unit_fix_v3:  Correspondances with Brian Walsh on 6/24/15 and
% 6/30/15, same as v2 except retains l/Lopt factor in second and third term
% of phi rate eq
clc, clear all, format long; 
tic

Ep = 5; % [J], pump energy
tp = 1e-3;  % [s], pump pulsewidth

wl_p = 790e-9;  % [m], pump wavelength
wl_l = 2.063e-6;  % [m], laser wavelength

Roc = 0.79;  % [dimensionless], output coupler reflectivity
R2 = 0.98;  % [dimensionless], passive cavity losses
Lc = 60e-2;  %  [m], cavity length
l = 4e-2;  % [m], rod length
w_rod = 0.2e-2;  % [m], rod radius
w_mode = 1.3e-3;  % [m], laser mode radius

n = 1.4641;  % [dimensionless], index of refraction of YLF at 2.063um
np = 1.3;  % [dimensionless], pump distribution factor (fudge factor, see paper)
eps = 2;  % [dimensionless], 2 for linear cavity, 1 for ring cavity
%Ts = 1;  % [dimensionless], laser rod surface transmission
passive_loss = 0;

Lopt = (n-1)*l + Lc;  % [m], cavity optical pathlength

Bspon = (pi*(w_rod^2))/(4*pi*(Lopt^2));
% Bspon = 1e-8;  % [dimensionless], spontaneous emission seed factor (fraction of photons spontaneously emitted within solid angle of cavity mirrors, value taken from "Computational model for operation of 2um co-doped Tm,Ho solid state lasers")
% Bspon = 0;

Ctm = 6/100;  % [dimensionless], Tm atomic concentration
Cho = 0.4/100;  % [dimensionless], Ho atomic concentration
rho_ylf = 3.99;  % [g/cm^3], YLF density
mY = 88.90585; % [amu], atomic mass of Yttrium
mLi = 6.941;  % [amu], atomic mass of Lithium
mF = 18.9984032;  % [amu], atomic mass of Fluorine
mp = 1.6726e-24;  % [g], proton mass
Ynum = 1;  % [dimensionless], number of Yttrium ions in YLF unit
Linum = 1;  % [dimensionless], number of Lithium ions in YLF unit
Fnum = 4;  % [dimensionless], number of Fluorine ions in YLF unit
n1_init = ((Ynum*rho_ylf)/((Ynum*mY + Linum*mLi + Fnum*mF)*mp))*Ctm*(1e6);  % [m^-3], Tm concentration density, assumed to all occupy ground manifold initially
n8_init = ((Ynum*rho_ylf)/((Ynum*mY + Linum*mLi + Fnum*mF)*mp))*Cho*(1e6);  % [m^-3], Ho concentration density, assumed to all occupy ground manifold initially

sigma_abs = (3.4e-21)*(1e-4);  % [m^2], pump absorption cross section (sigma polarization?? see paper)
sigma_em = (1.4e-19)*(1e-4);  % [m^2], stimulated emission cross section (sigma_em = sigma_se = sigma_e/f7,see paper) (pi polarization at 2.063um?? see paper)
% sigma_em = 0;

% energy transfer parameters (6%Tm, 0.4%Ho, from paper)
p28 = 1.68e-22;  % [m^3/s], Tm-Ho transfer
p71 = 1.28e-23;  % [m^3/s], Ho-Tm transfer
p41 = 3.13e-24;  % [m^3/s], Tm-Tm transfer
p22 = 3.48e-25;  % [m^3/s], Tm-Tm transfer
p38 = 7.64e-23;  % [m^3/s], Tm-Ho transfer
p61 = 4.09e-22;  % [m^3/s], Ho-Tm transfer
p27 = 1.47e-22;  % [m^3/s], Tm-Ho transfer
p51 = 1.16e-21;  % [m^3/s], Ho-Tm transfer
% p28 = 0;  % [m^3/s], Tm-Ho transfer
% p71 = 0;  % [m^3/s], Ho-Tm transfer
% p41 = 0;  % [m^3/s], Tm-Tm transfer
% p22 = 0;  % [m^3/s], Tm-Tm transfer
% p38 = 0;  % [m^3/s], Tm-Ho transfer
% p61 = 0;  % [m^3/s], Ho-Tm transfer
% p27 = 0;  % [m^3/s], Tm-Ho transfer
% p51 = 0;  % [m^3/s], Ho-Tm transfer

t2 = 15e-3;  % [s], Tm 3F4 lifetime
t3 = 1e-6;  % [s], Tm 3H5 lifetime
t4 = 2e-3;  % [s], Tm 3H4 lifetime
t5 = 20e-6;  % [s], Ho 5I5 lifetime
t6 = 2.2e-3;  % [s], Ho 5I6 lifetime
t7 = 16e-3;  % [s], Ho 5I7 lifetime

% branching ratios (see paper...no values given for B41 or B32, other values seem significantly different from those in literature from same authors)
% CLARIFICATION:  thest are really "branching probabilitie" that take into
% accountnon-radiative transitions to the next lower manifold, which for
% some cases can have much higher rates than the radiative transitions
% described by the explcit branching ratios
B41 = 0;
B32 = 1;
B43 = 1;
B56 = 1;
B67 = 1;

f7 = 0.0851;  % [dimensionless], Boltzmann fraction for upper laser level in YLF
f8 = 0.0258;  % [dimensionless], Boltzmann fraction for lower laser level in YLF

h = 6.626e-34;  % [J*s]
c = 2.998e8;  % [m/s]

Rp_level = (np*wl_p*Ep)/(h*c*pi*(w_rod^2)*l*tp);  % [m^-3*s^-1], pump rate
% tc = 1/((c/(eps*Lopt))*(log(Roc*R_passive)+2*(eps^2)*(1-Ts)));  % [s], cavity lifetime
tc = 1/((c/(eps*Lopt))*(-log(Roc)-log(R2)-log((1-passive_loss)^2)));  % [s], cavity lifetime (from "Computational model for operation of 2um co-doped Tm,Ho solid state lasers")
tr = (eps*Lopt)/c;  % [s], cavity roundtrip time

V_res=pi*((w_mode)^2)*Lopt;  % resonator mode volume  [m^3]
phi_init = 1/V_res;  % initial photon density = 1 photon per the resonator volume [m^-3]
%phi_init = 0;

%deltat = 0.1e-6;  %incremental time [s]
t_max = 10e-3;  %total cavity evolution time [s]

deltat = tp/1000;
tt = (0:deltat:t_max)';
Rp = (tt<=tp) * Rp_level;

%set up initial conditions for rate equations:
y_init = [n1_init,0,0,0,0,0,0,n8_init,phi_init,phi_init];
% y_init = [n1_init,0,0,0,0,0,0,n8_init];

%set up some tolerancing options for the ODE solver:
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'OutputFcn',@odeplot,'OutputSel',[1 2 4 7 8],'NonNegative',[1 2 3 4 5 6 7 8 9 10]);

%This is an array of parameters that is fed to the 'TmHoYLF_rate_equations.m' function and used in the following ODE solver:
p = [Rp_level,t2,t3,t4,t5,t6,t7,p28,p71,p41,p22,p38,p61,p27,p51,B41,B32,B43,B56,B67,f7,f8,c,Lopt,sigma_abs,sigma_em,tc,Bspon,Roc,tr,l,tp,n];

% Use ODE45 to solve coupled differential equations in 'TmHoYLF_rate_equations.m' as a function of time with the initial conditions and parameters defined above
figure(1);
[t,y,pump_rate] = ode15s('TmHoYLF_rate_equations_pulsed_pump_v2_test_unit_fix_v3',[0 t_max],y_init,options,p,tt,Rp);

% Here are the solutions to the coupled diff eqs solved by ODE45 [m^-3]
n1 = y(:,1);
n2 = y(:,2);
n3 = y(:,3);
n4 = y(:,4);
n5 = y(:,5);
n6 = y(:,6);
n7 = y(:,7);
n8 = y(:,8);
phi = y(:,9);
phi_out = y(:,10);

Eout_paper = ((h*c)/wl_l)*pi*(w_mode^2)*Lopt*phi_out(end);  %  [J], output energy using equation from paper

yTotal = y(:,1:8);
test = sum(yTotal');
yTm = y(:,1:4);
testTm = sum(yTm');
yHo = y(:,5:8);
testHo = sum(yHo');

figure(2); hold on
plot(t,test, 'g'); plot(t,testTm, 'r'); ...
plot(t,testHo, 'b'); plot(t,testTm+testHo,'--');
legend('total','Tm','Ho','Tm+Ho');
xlabel('time (s)');
ylabel('populations in manifolds (#)'); hold off


disp('                                                                ');
disp('**************************RESULTS*******************************');
disp('                                                                ');
disp(['Percentage of Tm atoms accounted for = ',num2str((sum(testTm)/(length(t)*n1_init))*100),'%']);
disp(['Percentage of Tm atoms accounted for = ',num2str((sum(testHo)/(length(t)*n8_init))*100),'%']);
disp(['Laser Output Pulse Energy = ',num2str(Eout_paper*1000),'mJ']);
disp('                                                                ');

figure(3);
subplot(4,1,1); plot(tt,Rp); legend('pump'); title(['v3 EndPump, R2=',num2str(R2)]);
subplot(4,1,2); plot(t,f7*n7,t,f8*n8); legend('f7n7','f8n8');
subplot(4,1,3); plot(t,f7*n7-f8*n8); legend('f7n7-f8n8');
subplot(4,1,4); plot(t,phi); legend('phi'); xlabel('time (s)')


toc

disp('                                                                ');